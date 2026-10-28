//! The `long` subcommand: trimming long reads (Oxford Nanopore, PacBio), following fastplong.
//!
//! Steps per read, in order (each can be disabled): fixed front/tail trimming, sliding-window
//! quality cutting, poly-X tail trimming, adapter trimming at the read start and end,
//! splitting at adapters inside the read (chimeras), removing low-quality stretches
//! (optional), and the filters, which are applied to every remaining piece. When a read is
//! split, its pieces are written as `<id>_1`, `<id>_2`, ... (numbered in read order before
//! filtering).
//!
//! Quality means (sliding windows, low-quality stretches, the mean quality filter and the
//! statistics) average error probabilities; fastplong averages phred values.

use std::cmp::Ordering;
use std::path::{Path, PathBuf};

use serde_json::{Value, json};

use crate::Result;
use crate::cli::{LongArgs, LowQualityMode};
use crate::fastq::{Chunk, Record};
use crate::filter::{LengthFilter, LongFate, LongFilter, LongQualityFilter, evaluate_long};
use crate::io::{self, OutputFiles};
use crate::log::Logger;
use crate::long_adapter::{self, LongAdapters};
use crate::long_detect::{Detection, Detector};
use crate::matcher::reverse_complement_into;
use crate::pipeline::{self, Sink};
use crate::platform::{Platform, parse_header};
use crate::polyx::poly_x_trim_len;
use crate::quality::{
    Averaging, QualityCut, Scale, Window, complement, low_quality_regions, mean_error_quality,
    trim_and_cut,
};
use crate::stats::{LongRunStats, ReadSetStats, Trimmed, group, percent, round};

/// Records per chunk at most.
const CHUNK_RECORDS: usize = 8192;
/// Chunks end once they hold this many bytes, so chunks of very long reads stay small.
const CHUNK_BYTES: usize = 8 << 20;

/// Output slots of a chunk result.
const OUT: usize = 0;
const FAILED: usize = 1;
const SLOTS: usize = 2;

/// Reads checked for `U` (RNA) and `T` bases before the run.
const BASE_CHECK_READS: usize = 100;

/// Where an adapter sequence came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Source {
    User,
    Detected,
    ReverseComplement,
}

impl Source {
    fn name(self) -> &'static str {
        match self {
            Source::User => "user",
            Source::Detected => "detected",
            Source::ReverseComplement => "reverse_complement_of_start",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum AdapterArg {
    Auto,
    None,
    Seq(Vec<u8>),
}

fn parse_adapter_arg(value: &str, what: &str) -> Result<AdapterArg> {
    match value.to_ascii_lowercase().as_str() {
        "auto" => return Ok(AdapterArg::Auto),
        "none" | "" => return Ok(AdapterArg::None),
        _ => {}
    }
    let seq: Vec<u8> = value
        .to_ascii_uppercase()
        .bytes()
        .map(|b| if b == b'U' { b'T' } else { b })
        .collect();
    if seq.len() <= 3 {
        return Err(format!("{what} must be longer than 3 bases"));
    }
    if let Some(&c) = seq.iter().find(|c| !b"ACGT".contains(c)) {
        return Err(format!(
            "{what} may only contain A, C, G, T or U (found '{}')",
            c as char
        ));
    }
    Ok(AdapterArg::Seq(seq))
}

/// Adapter sequences in use, for reporting.
struct AdapterInfo {
    start: Option<(Vec<u8>, Source)>,
    end: Option<(Vec<u8>, Source)>,
    fasta: usize,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Chimera {
    Off(&'static str),
    /// Split; keep pieces of at least this length.
    Split(usize),
    Discard,
}

struct LowQuality {
    mode: LowQualityMode,
    window: usize,
    threshold: u8,
}

/// Fully resolved trimming configuration.
struct Settings {
    front: usize,
    tail: usize,
    cut: QualityCut,
    poly_x: Option<usize>,
    adapters: Option<LongAdapters>,
    adapter_info: Option<AdapterInfo>,
    adapter_params: long_adapter::Params,
    chimera: Chimera,
    low_quality: Option<LowQuality>,
    filter: LongFilter,
    failed_output: bool,
    compression_level: u32,
    rna: bool,
}

fn validate_paths(args: &LongArgs) -> Result<()> {
    io::check_input(&args.input, "input")?;
    let mut inputs: Vec<&Path> = vec![&args.input];
    if let Some(fasta) = &args.adapter_fasta {
        io::check_input(fasta, "adapter FASTA")?;
        inputs.push(fasta);
    }
    let mut outputs: Vec<(&str, &Path)> = vec![("--output", &args.output)];
    if let Some(p) = &args.failed_output {
        outputs.push(("--failed-output", p));
    }
    if let Some(p) = &args.report {
        outputs.push(("--report", p));
    }
    io::check_outputs(&outputs, &inputs, args.force)
}

/// Looks at the first reads: header of the first read, and whether the reads use `U` (RNA).
fn inspect_start(path: &Path) -> Result<(Option<Vec<u8>>, bool)> {
    let mut input = io::open_fastq(path, "input")?;
    let Some(chunk) = input.reader.read_chunk(BASE_CHECK_READS, usize::MAX)? else {
        return Ok((None, false));
    };
    let (mut t, mut u) = (false, false);
    for rec in chunk.iter() {
        t |= rec.seq.contains(&b'T');
        u |= rec.seq.contains(&b'U');
    }
    if t && u {
        return Err(format!(
            "{}: the reads contain both U and T bases",
            path.display()
        ));
    }
    Ok((Some(chunk.get(0).name.to_vec()), u))
}

/// Runs adapter detection on the first reads of the input (fastplong's sample size).
fn detect_adapters(path: &Path, trim_tail: usize) -> Result<Detection> {
    let mut input = io::open_fastq(path, "input")?;
    let mut detector = Detector::new(trim_tail);
    while detector.wants_more() {
        let Some(chunk) = input.reader.read_chunk(4096, CHUNK_BYTES)? else {
            break;
        };
        for rec in chunk.iter() {
            if !detector.wants_more() {
                break;
            }
            detector.add(rec.seq);
        }
    }
    Ok(detector.finish())
}

fn resolve_adapters(args: &LongArgs, log: &Logger) -> Result<AdapterInfo> {
    let mut start = parse_adapter_arg(&args.start_adapter, "--start-adapter")?;
    let mut end = parse_adapter_arg(&args.end_adapter, "--end-adapter")?;
    let mut info = AdapterInfo {
        start: None,
        end: None,
        fasta: 0,
    };
    // as in fastplong: a given start adapter implies its reverse complement at the end
    if end == AdapterArg::Auto
        && let AdapterArg::Seq(s) = &start
    {
        let mut rc = Vec::new();
        reverse_complement_into(s, &mut rc);
        info.end = Some((rc, Source::ReverseComplement));
        end = AdapterArg::None;
    }
    if start == AdapterArg::Auto || end == AdapterArg::Auto {
        log.info("Detecting adapters ...");
        let d = detect_adapters(&args.input, args.trim_tail)?;
        for (which, arg, slot, flag, found) in [
            (
                "start",
                &mut start,
                &mut info.start,
                "--start-adapter",
                &d.start,
            ),
            ("end", &mut end, &mut info.end, "--end-adapter", &d.end),
        ] {
            if *arg != AdapterArg::Auto {
                continue;
            }
            *arg = AdapterArg::None;
            match &found.adapter {
                Some(seq) => {
                    log.info_with(
                        "adapter",
                        &format!("{which} adapter detected: {seq}"),
                        json!({ "read_end": which, "sequence": seq }),
                    );
                    let seq = seq
                        .bytes()
                        .map(|b| if b == b'U' { b'T' } else { b })
                        .collect();
                    *slot = Some((seq, Source::Detected));
                }
                None => {
                    let hint = found
                        .too_short
                        .as_ref()
                        .map(|s| {
                            format!(
                                " (frequent sequence {s} is too short; pass {flag} {s} if it is an adapter)"
                            )
                        })
                        .unwrap_or_default();
                    log.info_with(
                        "adapter",
                        &format!("{which}: no adapter detected{hint}"),
                        json!({ "read_end": which, "sequence": null, "too_short_candidate": found.too_short }),
                    );
                }
            }
        }
        if d.reads < crate::long_detect::MIN_READS {
            log.warn(&format!(
                "too few reads for adapter detection ({} of {} needed)",
                d.reads,
                crate::long_detect::MIN_READS
            ));
        }
    }
    if let AdapterArg::Seq(s) = start {
        info.start = Some((s, Source::User));
    }
    if let AdapterArg::Seq(s) = end {
        info.end = Some((s, Source::User));
    }
    Ok(info)
}

fn window(
    enabled: bool,
    size: Option<u16>,
    quality: Option<u8>,
    args: &LongArgs,
) -> Option<Window> {
    enabled.then(|| Window {
        size: size.unwrap_or(args.cut_window_size) as usize,
        quality: quality.unwrap_or(args.cut_mean_quality),
    })
}

fn resolve_settings(args: &LongArgs, log: &Logger) -> Result<Settings> {
    let (first_name, rna) = inspect_start(&args.input)?;
    if let Some(h) = first_name.as_deref().map(parse_header)
        && matches!(
            h.platform,
            Platform::Illumina | Platform::ElementAviti | Platform::Mgi
        )
    {
        log.warn(&format!(
            "the input looks like {} short reads; 'long' is meant for long reads ('simpletrim short' handles short reads)",
            h.platform.label()
        ));
    }

    let cut = QualityCut {
        front: window(
            args.cut_front,
            args.cut_front_window_size,
            args.cut_front_mean_quality,
            args,
        ),
        tail: window(
            args.cut_tail,
            args.cut_tail_window_size,
            args.cut_tail_mean_quality,
            args,
        ),
        right: None,
    };

    let adapter_params = long_adapter::Params {
        rate: args.distance_threshold,
        min_overlap: args.adapter_min_overlap as usize,
        window: args.adapter_window as usize,
        extension: args.trimming_extension as usize,
        internal_rate: args.chimera_distance_threshold,
    };
    let (adapters, adapter_info) = if args.disable_adapter_trimming {
        (None, None)
    } else {
        let mut info = resolve_adapters(args, log)?;
        let fasta = match &args.adapter_fasta {
            Some(p) => io::load_adapter_fasta(p, log)?,
            None => Vec::new(),
        };
        info.fasta = fasta.len();
        let adapters =
            (info.start.is_some() || info.end.is_some() || !fasta.is_empty()).then(|| {
                LongAdapters::new(
                    info.start.as_ref().map(|(s, _)| s.as_slice()),
                    info.end.as_ref().map(|(s, _)| s.as_slice()),
                    &fasta,
                    adapter_params,
                )
            });
        (adapters, Some(info))
    };

    let chimera = match &adapters {
        None => Chimera::Off("no adapter sequence"),
        Some(_) if args.disable_chimera_splitting => Chimera::Off("disabled"),
        Some(a) if !a.has_internal() => {
            Chimera::Off("no adapter long enough to search inside reads")
        }
        Some(_) if args.discard_chimeric_reads => Chimera::Discard,
        Some(_) => Chimera::Split(args.chimera_min_segment_length as usize),
    };

    let filter = LongFilter {
        quality: (!args.disable_quality_filter).then_some(LongQualityFilter {
            qualified_quality: args.qualified_quality,
            unqualified_percent_limit: args.unqualified_percent_limit,
            n_percent_limit: args.n_percent_limit,
            n_base_limit: args.n_base_limit,
            average_quality: args.average_quality,
        }),
        length: (!args.disable_length_filter).then_some(LengthFilter {
            min_len: args.min_length,
            max_len: args.max_length,
        }),
        gc: (args.min_gc.is_some() || args.max_gc.is_some())
            .then(|| (args.min_gc.unwrap_or(0.0), args.max_gc.unwrap_or(100.0))),
        complexity: args
            .low_complexity_filter
            .then_some(args.complexity_threshold as f64 / 100.0),
    };
    if filter.length.is_some() && args.max_length > 0 && args.max_length < args.min_length {
        return Err("--max-length is smaller than --min-length".into());
    }
    if let Some((min, max)) = filter.gc
        && min > max
    {
        return Err("--min-gc is larger than --max-gc".into());
    }

    Ok(Settings {
        front: args.trim_front,
        tail: args.trim_tail,
        cut,
        poly_x: args.trim_poly_x.then_some(args.poly_x_min_len as usize),
        adapters,
        adapter_info,
        adapter_params,
        chimera,
        low_quality: args.low_quality_mode.map(|mode| LowQuality {
            mode,
            window: args.low_quality_window as usize,
            threshold: args.low_quality_threshold,
        }),
        filter,
        failed_output: args.failed_output.is_some(),
        compression_level: args.compression_level,
        rna,
    })
}

/// The longest piece; on ties the one with the better mean quality, then the first.
fn longest(pieces: &[(usize, usize)], qual: &[u8]) -> Option<(usize, usize)> {
    let mean = |(a, b): (usize, usize)| mean_error_quality(&qual[a..b]).unwrap_or(0.0);
    pieces.iter().copied().max_by(|&x, &y| {
        (x.1 - x.0)
            .cmp(&(y.1 - y.0))
            .then_with(|| mean(x).partial_cmp(&mean(y)).unwrap_or(Ordering::Equal))
            .then_with(|| y.0.cmp(&x.0))
    })
}

/// Shifts ranges relative to `offset` into read coordinates.
fn shifted(ranges: Vec<(usize, usize)>, offset: usize) -> Vec<(usize, usize)> {
    ranges
        .into_iter()
        .map(|(a, b)| (a + offset, b + offset))
        .collect()
}

fn process_read(
    s: &Settings,
    rec: &Record,
    st: &mut LongRunStats,
    out: &mut [Vec<u8>; SLOTS],
    buf: &mut Vec<u8>,
) {
    let errp = Scale::get(Averaging::ErrorProbability);
    st.input.add(rec.seq, rec.qual);
    let len = rec.seq.len();
    let kept = trim_and_cut(rec.seq, rec.qual, s.front, s.tail, &s.cut, errp);
    let fixed = (s.front + s.tail).min(len);
    st.trim.fixed.add(fixed);
    st.trim
        .quality_window
        .add((len - fixed).saturating_sub(kept.map_or(0, |(a, b)| b - a)));

    let mut pieces: Vec<(usize, usize)> = Vec::new();
    // what happened to a read that leaves nothing
    let mut empty_fate = LongFate::TooShort;
    if let Some((mut a, mut b)) = kept {
        if let Some(min_len) = s.poly_x
            && let Some(new_len) = poly_x_trim_len(&rec.seq[a..b], min_len)
        {
            st.trim.poly_x.add(b - a - new_len);
            b = a + new_len;
        }
        if let Some(ad) = &s.adapters {
            let cut = ad.trim_start(&rec.seq[a..b]);
            st.trim.adapter_start.add(cut);
            a += cut;
            let cut = ad.trim_end(&rec.seq[a..b], buf);
            st.trim.adapter_end.add(cut);
            b -= cut;
        }
        if b > a {
            pieces.push((a, b));
        }

        if let (Some(ad), Some(&(a, b))) = (&s.adapters, pieces.first())
            && !matches!(s.chimera, Chimera::Off(_))
        {
            let regions = ad.internal_regions(&rec.seq[a..b]);
            if !regions.is_empty() {
                st.chimeras.reads += 1;
                if let Chimera::Split(min_len) = s.chimera {
                    let parts: Vec<(usize, usize)> = shifted(complement(b - a, &regions), a)
                        .into_iter()
                        .filter(|&(x, y)| y - x >= min_len)
                        .collect();
                    let kept_bases: usize = parts.iter().map(|&(x, y)| y - x).sum();
                    st.chimeras.bases_removed += (b - a - kept_bases) as u64;
                    st.chimeras.pieces += parts.len() as u64;
                    pieces = parts;
                } else {
                    st.chimeras.reads_discarded += 1;
                    pieces.clear();
                }
                if pieces.is_empty() {
                    empty_fate = LongFate::Chimeric;
                }
            }
        }

        if let Some(lq) = &s.low_quality
            && !pieces.is_empty()
        {
            let mut next = Vec::with_capacity(pieces.len());
            let mut affected = false;
            let mut removed = 0;
            for &(a, b) in &pieces {
                let regions = low_quality_regions(&rec.qual[a..b], lq.window, lq.threshold, errp);
                if regions.is_empty() {
                    next.push((a, b));
                    continue;
                }
                affected = true;
                let good = shifted(complement(b - a, &regions), a);
                let chosen = match lq.mode {
                    LowQualityMode::Split => good,
                    LowQualityMode::Longest => longest(&good, rec.qual).into_iter().collect(),
                };
                removed += (b - a) - chosen.iter().map(|&(x, y)| y - x).sum::<usize>();
                next.extend(chosen);
            }
            if affected {
                st.low_quality.reads += 1;
                st.low_quality.bases_removed += removed as u64;
                st.low_quality.pieces += next.len() as u64;
                if next.is_empty() {
                    empty_fate = LongFate::LowQuality;
                }
            }
            pieces = next;
        }
    }

    let numbered = pieces.len() > 1;
    let mut written = 0;
    for (k, &(a, b)) in pieces.iter().enumerate() {
        let fate = evaluate_long(&rec.seq[a..b], &rec.qual[a..b], &s.filter);
        st.fates.add(fate);
        let piece = numbered.then_some(k + 1);
        if fate == LongFate::Pass {
            rec.write_piece(a, b, piece, None, &mut out[OUT]);
            st.output.add(&rec.seq[a..b], &rec.qual[a..b]);
            written += 1;
        } else if s.failed_output {
            rec.write_piece(a, b, piece, Some(fate.tag()), &mut out[FAILED]);
        }
    }
    if pieces.is_empty() {
        st.fates.add(empty_fate);
        if s.failed_output {
            rec.write_piece(0, len, None, Some(empty_fate.tag()), &mut out[FAILED]);
        }
    }
    st.reads_with_output += (written > 0) as u64;
    st.reads_split += (written > 1) as u64;
}

fn process_chunk(s: &Settings, chunk: &Chunk) -> (Vec<Vec<u8>>, LongRunStats) {
    let mut st = LongRunStats::default();
    let mut out: [Vec<u8>; SLOTS] = Default::default();
    let mut buf = Vec::new();
    for rec in chunk.iter() {
        process_read(s, &rec, &mut st, &mut out, &mut buf);
    }
    (out.into(), st)
}

fn describe_settings(s: &Settings) -> (Vec<String>, Value) {
    let mut lines = Vec::new();
    lines.push(if s.front > 0 || s.tail > 0 {
        format!("fixed trimming: front {} tail {}", s.front, s.tail)
    } else {
        "fixed trimming: off".into()
    });
    let win = |name: &str, w: Option<Window>| {
        w.map(|w| format!("{name} (window {}, Q{})", w.size, w.quality))
    };
    let cuts: Vec<String> = [win("front", s.cut.front), win("tail", s.cut.tail)]
        .into_iter()
        .flatten()
        .collect();
    lines.push(if cuts.is_empty() {
        "quality window: off".into()
    } else {
        format!(
            "quality window (error-probability mean): {}",
            cuts.join(", ")
        )
    });
    lines.push(match s.poly_x {
        Some(n) => format!("poly-X trimming: on, min length {n}"),
        None => "poly-X trimming: off".into(),
    });
    let p = &s.adapter_params;
    let mut adapter_json = Value::Null;
    match &s.adapter_info {
        None => lines.push("adapter trimming: off".into()),
        Some(info) => {
            let mut parts = Vec::new();
            for (name, a) in [("start", &info.start), ("end", &info.end)] {
                if let Some((seq, src)) = a {
                    parts.push(format!(
                        "{name} {} ({})",
                        String::from_utf8_lossy(seq),
                        src.name()
                    ));
                }
            }
            if info.fasta > 0 {
                parts.push(format!("{} FASTA sequences", info.fasta));
            }
            if parts.is_empty() {
                lines.push("adapter trimming: on, but no adapter sequence".into());
            } else {
                lines.push(format!(
                    "adapter trimming: {}; up to {}% edits, min overlap {}, window {} bp, +{} bp",
                    parts.join(" + "),
                    round(p.rate * 100.0, 1),
                    p.min_overlap,
                    p.window,
                    p.extension
                ));
            }
            let seq_json = |a: &Option<(Vec<u8>, Source)>| {
                json!({
                    "sequence": a.as_ref().map(|(s, _)| String::from_utf8_lossy(s).into_owned()),
                    "source": a.as_ref().map(|(_, src)| src.name()),
                })
            };
            adapter_json = json!({
                "start": seq_json(&info.start),
                "end": seq_json(&info.end),
                "fasta_sequences": info.fasta,
                "distance_threshold": p.rate,
                "min_overlap": p.min_overlap,
                "window": p.window,
                "trimming_extension": p.extension,
            });
        }
    }
    let chimera_json = match s.chimera {
        Chimera::Off(reason) => {
            lines.push(format!("chimeric reads: not checked ({reason})"));
            json!({ "mode": "off", "reason": reason })
        }
        Chimera::Split(min_len) => {
            lines.push(format!(
                "chimeric reads: split at internal adapters (up to {}% edits), pieces >= {min_len} bp",
                round(p.internal_rate * 100.0, 1)
            ));
            json!({ "mode": "split", "min_segment_length": min_len, "distance_threshold": p.internal_rate })
        }
        Chimera::Discard => {
            lines.push(format!(
                "chimeric reads: discarded (internal adapters with up to {}% edits)",
                round(p.internal_rate * 100.0, 1)
            ));
            json!({ "mode": "discard", "distance_threshold": p.internal_rate })
        }
    };
    lines.push(match &s.low_quality {
        Some(lq) => format!(
            "low-quality stretches: {} (window {}, mean below Q{})",
            lq.mode.name(),
            lq.window,
            lq.threshold
        ),
        None => "low-quality stretches: kept".into(),
    });
    let f = &s.filter;
    lines.push(match &f.quality {
        Some(q) => format!(
            "quality filter: <= {}% bases below Q{}, <= {}% N{}{}",
            q.unqualified_percent_limit,
            q.qualified_quality,
            q.n_percent_limit,
            q.n_base_limit
                .map(|n| format!(", <= {n} N"))
                .unwrap_or_default(),
            if q.average_quality > 0 {
                format!(", mean >= Q{}", q.average_quality)
            } else {
                String::new()
            }
        ),
        None => "quality filter: off".into(),
    });
    lines.push(match &f.length {
        Some(l) if l.max_len > 0 => format!("length filter: {}-{} bp", l.min_len, l.max_len),
        Some(l) => format!("length filter: >= {} bp", l.min_len),
        None => "length filter: off".into(),
    });
    if let Some((min, max)) = f.gc {
        lines.push(format!("GC filter: {min}-{max}%"));
    }
    if let Some(c) = f.complexity {
        lines.push(format!("complexity filter: >= {}%", round(c * 100.0, 1)));
    }
    if s.rna {
        lines.push("reads contain U (RNA)".into());
    }
    let json = json!({
        "trim_front": s.front,
        "trim_tail": s.tail,
        "cut_front": s.cut.front.map(|w| json!({"window": w.size, "quality": w.quality})),
        "cut_tail": s.cut.tail.map(|w| json!({"window": w.size, "quality": w.quality})),
        "quality_averaging": "error_probability",
        "poly_x": { "enabled": s.poly_x.is_some(), "min_len": s.poly_x },
        "adapter": adapter_json,
        "chimeras": chimera_json,
        "low_quality": s.low_quality.as_ref().map(|lq| json!({
            "mode": lq.mode.name(),
            "window": lq.window,
            "threshold": lq.threshold,
        })),
        "quality_filter": f.quality.map(|q| json!({
            "qualified_quality": q.qualified_quality,
            "unqualified_percent_limit": q.unqualified_percent_limit,
            "n_percent_limit": q.n_percent_limit,
            "n_base_limit": q.n_base_limit,
            "average_quality": q.average_quality,
        })),
        "length_filter": f.length.map(|l| json!({"min_len": l.min_len, "max_len": l.max_len})),
        "gc_filter": f.gc.map(|(min, max)| json!({"min": min, "max": max})),
        "complexity_filter": f.complexity,
        "rna": s.rna,
        "compression_level": s.compression_level,
    });
    (lines, json)
}

fn set_line(label: &str, r: &ReadSetStats) -> String {
    let fmt_opt = |v: Option<usize>| v.map_or("-".into(), |n| group(n as u64));
    format!(
        "{label} {} reads, {} bases, mean length {:.1}, N50 {}, mean quality {}",
        group(r.base.reads),
        group(r.base.bases),
        r.base.mean_length(),
        fmt_opt(r.lengths.nx(0.5)),
        r.mean_quality().map_or("-".into(), |q| format!("Q{q:.1}"))
    )
}

fn summary(s: &Settings, args: &LongArgs, st: &LongRunStats, elapsed: f64) -> (String, Value) {
    let mut t = String::new();
    let mut line = |l: String| {
        t.push_str(&l);
        t.push('\n');
    };
    let reads = st.input.base.reads;
    line(format!("Done in {elapsed:.1} s."));
    line(set_line("Input: ", &st.input));
    line(format!(
        "Output: {} records from {} reads ({} of reads, {} of bases)",
        group(st.fates.passed),
        group(st.reads_with_output),
        percent(st.reads_with_output, reads),
        percent(st.output.base.bases, st.input.base.bases)
    ));
    line(set_line("       ", &st.output));
    if st.reads_split > 0 {
        line(format!(
            "  reads written as several records: {}",
            group(st.reads_split)
        ));
    }
    line("Trimmed (reads / bases):".into());
    let tr = &st.trim;
    let cell = |x: &Trimmed| format!("{:>12} / {}", group(x.reads), group(x.bases));
    for (name, x) in [
        ("fixed", &tr.fixed),
        ("quality window", &tr.quality_window),
        ("poly-X", &tr.poly_x),
        ("adapter (start)", &tr.adapter_start),
        ("adapter (end)", &tr.adapter_end),
    ] {
        line(format!("  {name:<16} {}", cell(x)));
    }
    let c = &st.chimeras;
    match s.chimera {
        Chimera::Off(_) => {}
        Chimera::Split(_) => line(format!(
            "Chimeric reads: {} with internal adapters, split into {} pieces ({} bases removed)",
            group(c.reads),
            group(c.pieces),
            group(c.bases_removed)
        )),
        Chimera::Discard => line(format!(
            "Chimeric reads: {} with internal adapters, discarded",
            group(c.reads)
        )),
    }
    if s.low_quality.is_some() {
        let l = &st.low_quality;
        line(format!(
            "Low-quality stretches: {} reads, {} bases removed, {} pieces kept",
            group(l.reads),
            group(l.bases_removed),
            group(l.pieces)
        ));
    }
    line("Record fates (reads or pieces):".into());
    let f = &st.fates;
    for (name, n) in [
        ("passed", f.passed),
        ("too short", f.too_short),
        ("too long", f.too_long),
        ("low quality", f.low_quality),
        ("too many N", f.too_many_n),
        ("GC content", f.gc_content),
        ("low complexity", f.low_complexity),
        ("chimeric", f.chimeric),
    ] {
        line(format!("  {name:<16} {:>14}", group(n)));
    }

    let path = |p: &Option<PathBuf>| p.as_ref().map(|p| p.display().to_string());
    let json = json!({
        "mode": "long",
        "elapsed_s": round(elapsed, 3),
        "inputs": { "input": args.input.display().to_string() },
        "outputs": {
            "output": args.output.display().to_string(),
            "failed_output": path(&args.failed_output),
        },
        "reads_in": reads,
        "records_written": st.fates.passed,
        "reads_with_output": st.reads_with_output,
        "reads_split": st.reads_split,
        "before": st.input.to_json(),
        "after": st.output.to_json(),
        "trimming": serde_json::to_value(st.trim).expect("serializable"),
        "chimeras": if matches!(s.chimera, Chimera::Off(_)) { Value::Null } else { serde_json::to_value(st.chimeras).expect("serializable") },
        "low_quality": if s.low_quality.is_some() { serde_json::to_value(st.low_quality).expect("serializable") } else { Value::Null },
        "filtering": serde_json::to_value(st.fates).expect("serializable"),
    });
    (t, json)
}

/// Runs the `long` subcommand.
pub fn run(args: &LongArgs, threads: usize, log: &Logger) -> Result<()> {
    validate_paths(args)?;
    let settings = resolve_settings(args, log)?;
    let (lines, settings_json) = describe_settings(&settings);
    let mut header = format!(
        "simpletrim {} long | {threads} worker threads\n",
        env!("CARGO_PKG_VERSION")
    );
    for l in &lines {
        header.push_str("  ");
        header.push_str(l);
        header.push('\n');
    }
    log.info_with(
        "start",
        header.trim_end(),
        json!({
            "version": env!("CARGO_PKG_VERSION"),
            "mode": "long",
            "threads": threads,
            "settings": settings_json,
        }),
    );

    let mut files = OutputFiles::new();
    let mut sinks: Vec<Option<Sink>> = (0..SLOTS).map(|_| None).collect();
    sinks[OUT] = Some(Sink::create(&mut files, &args.output, args.force)?);
    if let Some(p) = &args.failed_output {
        sinks[FAILED] = Some(Sink::create(&mut files, p, args.force)?);
    }
    let options = pipeline::Options {
        threads,
        max_reads: args.max_reads,
        chunk_records: CHUNK_RECORDS,
        chunk_bytes: CHUNK_BYTES,
        progress_interval: args.progress_interval,
        unit: "reads",
        compression_level: settings.compression_level,
        log,
    };
    let input = io::open_fastq(&args.input, "input")?;
    let stats = pipeline::run(&options, input, None, sinks, |chunk, _| {
        Ok(process_chunk(&settings, chunk))
    })?;
    let elapsed = log.elapsed_secs();
    let (text, mut json) = summary(&settings, args, &stats, elapsed);
    json["settings"] = settings_json;
    if let Some(report) = &args.report {
        let pretty = serde_json::to_string_pretty(&json).expect("serializable");
        io::write_text_file(report, &(pretty + "\n"), args.force)?;
    }
    files.commit();
    log.result("summary", &text, json);
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::cli::{Cli, Command};
    use crate::fastq::FastqReader;
    use clap::Parser;

    fn args(extra: &[&str]) -> LongArgs {
        let mut argv = vec!["simpletrim", "long", "-i", "in.fq", "-o", "out.fq"];
        argv.extend_from_slice(extra);
        match Cli::try_parse_from(argv).unwrap().command {
            Command::Long(a) => *a,
            _ => unreachable!(),
        }
    }

    /// Settings without file access: adapters given explicitly.
    fn settings(a: &LongArgs, start: &[u8], end: &[u8]) -> Settings {
        let params = long_adapter::Params {
            rate: a.distance_threshold,
            min_overlap: a.adapter_min_overlap as usize,
            window: a.adapter_window as usize,
            extension: a.trimming_extension as usize,
            internal_rate: a.chimera_distance_threshold,
        };
        let adapters = LongAdapters::new(Some(start), Some(end), &[], params);
        Settings {
            front: a.trim_front,
            tail: a.trim_tail,
            cut: QualityCut::default(),
            poly_x: None,
            adapters: Some(adapters),
            adapter_info: None,
            adapter_params: params,
            chimera: if a.discard_chimeric_reads {
                Chimera::Discard
            } else {
                Chimera::Split(a.chimera_min_segment_length as usize)
            },
            low_quality: a.low_quality_mode.map(|mode| LowQuality {
                mode,
                window: a.low_quality_window as usize,
                threshold: a.low_quality_threshold,
            }),
            filter: LongFilter {
                quality: None,
                length: Some(LengthFilter {
                    min_len: a.min_length,
                    max_len: 0,
                }),
                gc: None,
                complexity: None,
            },
            failed_output: true,
            compression_level: 1,
            rna: false,
        }
    }

    fn run_reads(s: &Settings, fastq: &str) -> (String, String, LongRunStats) {
        let mut r = FastqReader::new(
            Box::new(std::io::Cursor::new(fastq.as_bytes().to_vec())),
            "test".into(),
        );
        let chunk = r.read_chunk(1000, usize::MAX).unwrap().unwrap();
        let (out, st) = process_chunk(s, &chunk);
        (
            String::from_utf8(out[OUT].clone()).unwrap(),
            String::from_utf8(out[FAILED].clone()).unwrap(),
            st,
        )
    }

    const START: &[u8] = b"ACTTCGTTCAGTTACGTATTGCTGGATCCAGTA";
    const END: &[u8] = b"GCAATACGTAACTGAACGAAGTACAGGTTCA";

    fn insert(seed: u64, len: usize) -> Vec<u8> {
        crate::synth::Rng::new(seed).seq(len)
    }

    fn fq(name: &str, seq: &[u8], qual: &[u8]) -> String {
        format!(
            "@{name}\n{}\n+\n{}\n",
            String::from_utf8_lossy(seq),
            String::from_utf8_lossy(qual)
        )
    }

    #[test]
    fn trims_ends_and_splits_chimeras() {
        let a = args(&[]);
        let s = settings(&a, START, END);
        let (i1, i2) = (insert(1, 300), insert(2, 200));
        let mut seq = START.to_vec();
        seq.extend_from_slice(&i1);
        seq.extend_from_slice(END);
        seq.extend_from_slice(START);
        seq.extend_from_slice(&i2);
        seq.extend_from_slice(&END[..20]);
        let qual = vec![b'5'; seq.len()];
        let (out, failed, st) = run_reads(&s, &fq("r1 comment", &seq, &qual));
        let recs: Vec<&str> = out.lines().collect();
        assert_eq!(recs.len(), 8, "{out}");
        assert_eq!(recs[0], "@r1_1 comment");
        assert_eq!(recs[4], "@r1_2 comment");
        // the pieces are the inserts without the 10 bp extensions
        assert_eq!(recs[1].as_bytes(), &i1[10..i1.len() - 10]);
        assert_eq!(recs[5].as_bytes(), &i2[10..i2.len() - 10]);
        assert!(failed.is_empty());
        assert_eq!((st.chimeras.reads, st.chimeras.pieces), (1, 2));
        assert_eq!((st.reads_with_output, st.reads_split), (1, 1));
        assert_eq!(st.trim.adapter_start.bases, START.len() as u64 + 10);
    }

    #[test]
    fn discard_mode_and_failed_output() {
        let a = args(&["--discard-chimeric-reads"]);
        let s = settings(&a, START, END);
        let mut chimera = insert(3, 200);
        chimera.extend_from_slice(START);
        chimera.extend(insert(4, 200));
        let short = insert(5, 15);
        let mut text = fq("chim", &chimera, &vec![b'5'; chimera.len()]);
        text.push_str(&fq("short", &short, &vec![b'5'; short.len()]));
        let (out, failed, st) = run_reads(&s, &text);
        assert!(out.is_empty());
        let names: Vec<&str> = failed.lines().step_by(4).collect();
        assert_eq!(
            names,
            ["@chim failed_chimeric_read", "@short failed_too_short"]
        );
        assert_eq!(st.fates.chimeric, 1);
        assert_eq!(st.fates.too_short, 1);
    }

    #[test]
    fn low_quality_modes() {
        let seq = insert(6, 1000);
        let mut qual = vec![b'5'; 1000];
        for q in &mut qual[300..420] {
            *q = b'#'; // Q2
        }
        for q in &mut qual[800..860] {
            *q = b'#';
        }
        let text = fq("r", &seq, &qual);
        let a = args(&["--low-quality-mode", "split"]);
        let (out, _, st) = run_reads(&settings(&a, START, END), &text);
        let lens: Vec<usize> = out.lines().skip(1).step_by(4).map(str::len).collect();
        assert_eq!(lens, [300, 380, 140]);
        assert_eq!(st.low_quality.bases_removed, 180);
        let a = args(&["--low-quality-mode", "longest"]);
        let (out, _, _) = run_reads(&settings(&a, START, END), &text);
        let recs: Vec<&str> = out.lines().collect();
        assert_eq!(recs.len(), 4);
        assert_eq!(recs[0], "@r");
        assert_eq!(recs[1].as_bytes(), &seq[420..800]);
    }

    #[test]
    fn adapter_arguments() {
        assert_eq!(parse_adapter_arg("auto", "x"), Ok(AdapterArg::Auto));
        assert_eq!(parse_adapter_arg("None", "x"), Ok(AdapterArg::None));
        assert_eq!(
            parse_adapter_arg("acgu", "x"),
            Ok(AdapterArg::Seq(b"ACGT".to_vec()))
        );
        assert!(parse_adapter_arg("ACG", "x").is_err());
        assert!(parse_adapter_arg("ACGN", "x").is_err());
    }
}
