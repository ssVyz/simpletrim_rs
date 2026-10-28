//! The `short` subcommand: trimming Illumina-style single-end and paired-end reads.
//!
//! Steps per read, in order (each can be disabled): fixed front/tail trimming, sliding-window
//! quality cutting, poly-G tail trimming, adapter trimming (read overlap for pairs, then
//! adapter sequences), adapter dimer check, quality and length filters.
//!
//! Data flow: one reader thread per input parses fixed-size chunks; chunks are processed and
//! (for .gz outputs) compressed in parallel on the worker pool; a writer thread writes the
//! results in input order.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, SyncSender, channel, sync_channel};
use std::time::{Duration, Instant};

use serde_json::{Value, json};

use crate::Result;
use crate::adapter::{MatchParams, trim_multi, trim_single};
use crate::cli::{OrphanMode, ShortArgs};
use crate::detect::{self, Detected};
use crate::fastq::{Chunk, Record, feed_chunks, mate_names_match};
use crate::filter::{self, Fate, LengthFilter, QualityFilter};
use crate::io::{self, OutputFiles, gzip_member, is_gz};
use crate::known_adapters;
use crate::log::Logger;
use crate::overlap::{self, OverlapParams};
use crate::platform::{Platform, parse_header};
use crate::polyg::poly_g_trim_len;
use crate::quality::{QualityCut, Window, trim_and_cut};
use crate::stats::{RunStats, TrimStats, group, percent};

/// Records per chunk; both mates are read in chunks of the same size so pairs stay aligned.
const CHUNK_RECORDS: usize = 8192;

/// Output slots of a chunk result.
const OUT1: usize = 0;
const OUT2: usize = 1;
const ORPHAN1: usize = 2;
const ORPHAN2: usize = 3;
const SLOTS: usize = 4;

/// Where an adapter sequence came from.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Source {
    User,
    Preset,
    Detected,
    SameAsR1,
}

impl Source {
    fn name(self) -> &'static str {
        match self {
            Source::User => "user",
            Source::Preset => "preset",
            Source::Detected => "detected",
            Source::SameAsR1 => "same_as_r1",
        }
    }
}

#[derive(Clone, Debug, PartialEq, Eq)]
enum AdapterArg {
    Auto,
    None,
    Seq(Vec<u8>, Source),
}

fn parse_adapter_arg(value: &str, what: &str) -> Result<AdapterArg> {
    match value.to_ascii_lowercase().as_str() {
        "auto" => return Ok(AdapterArg::Auto),
        "none" | "" => return Ok(AdapterArg::None),
        _ => {}
    }
    let seq = value.to_ascii_uppercase().into_bytes();
    if seq.len() <= 3 {
        return Err(format!("{what} must be longer than 3 bases"));
    }
    if let Some(&c) = seq.iter().find(|c| !b"ACGT".contains(c)) {
        return Err(format!(
            "{what} may only contain A, C, G, T (found '{}')",
            c as char
        ));
    }
    Ok(AdapterArg::Seq(seq, Source::User))
}

/// Adapter trimming configuration.
pub struct AdapterSettings {
    seq: [Option<Vec<u8>>; 2],
    source: [Option<Source>; 2],
    names: [Option<&'static str>; 2],
    fasta: Vec<Vec<u8>>,
    matching: MatchParams,
    overlap: OverlapParams,
}

/// Fully resolved trimming configuration.
pub struct Settings {
    paired: bool,
    front: [usize; 2],
    tail: [usize; 2],
    cut: QualityCut,
    poly_g: Option<usize>,
    poly_g_reason: &'static str,
    adapter: Option<AdapterSettings>,
    dimer_max_len: Option<usize>,
    quality_filter: Option<QualityFilter>,
    length_filter: Option<LengthFilter>,
    orphan_mode: OrphanMode,
    check_names: bool,
    /// Output slot that receives read 2 orphans.
    orphan2_slot: usize,
    /// Which output slots are gzip-compressed.
    gz: [bool; SLOTS],
    compression_level: u32,
}

fn window(
    enabled: bool,
    size: Option<u16>,
    quality: Option<u8>,
    args: &ShortArgs,
) -> Option<Window> {
    enabled.then(|| Window {
        size: size.unwrap_or(args.cut_window_size) as usize,
        quality: quality.unwrap_or(args.cut_mean_quality),
    })
}

/// Reads a FASTA file of adapters (sequences of at least 6 bases, de-duplicated, as fastp).
fn load_adapter_fasta(path: &Path, log: &Logger) -> Result<Vec<Vec<u8>>> {
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("cannot read adapter FASTA {}: {e}", path.display()))?;
    let mut records: Vec<(String, Vec<u8>)> = Vec::new();
    for line in text.lines().map(str::trim) {
        if let Some(name) = line.strip_prefix('>') {
            records.push((name.trim().to_string(), Vec::new()));
        } else if !line.is_empty() {
            let Some((name, seq)) = records.last_mut() else {
                return Err(format!(
                    "adapter FASTA {}: sequence before the first '>' header",
                    path.display()
                ));
            };
            let upper = line.to_ascii_uppercase();
            if let Some(c) = upper.chars().find(|c| !"ACGTN".contains(*c)) {
                return Err(format!(
                    "adapter FASTA {}: invalid base '{c}' in '{name}'",
                    path.display()
                ));
            }
            seq.extend_from_slice(upper.as_bytes());
        }
    }
    let mut adapters: Vec<Vec<u8>> = Vec::new();
    for (name, seq) in records {
        if seq.len() < 6 {
            log.warn(&format!(
                "skipping adapter '{name}' from {}: shorter than 6 bases",
                path.display()
            ));
        } else if !adapters.contains(&seq) {
            adapters.push(seq);
        }
    }
    if adapters.is_empty() {
        return Err(format!(
            "adapter FASTA {} has no usable sequences",
            path.display()
        ));
    }
    Ok(adapters)
}

fn peek_first_name(path: &Path, label: &str) -> Result<Option<Vec<u8>>> {
    let mut input = io::open_fastq(path, label)?;
    Ok(input
        .reader
        .read_chunk(1, usize::MAX)?
        .map(|c| c.get(0).name.to_vec()))
}

fn validate_paths(args: &ShortArgs) -> Result<()> {
    io::check_input(&args.r1, "R1 input")?;
    let mut inputs: Vec<&Path> = vec![&args.r1];
    if let Some(r2) = &args.r2 {
        io::check_input(r2, "R2 input")?;
        if io::same_path(&args.r1, r2) {
            return Err("--r1 and --r2 point to the same file".into());
        }
        inputs.push(r2);
    }
    if let Some(fasta) = &args.adapter_fasta {
        io::check_input(fasta, "adapter FASTA")?;
        inputs.push(fasta);
    }
    let orphan_files = args.orphan1.is_some() || args.orphan2.is_some();
    match (args.orphan_mode, args.r2.is_some()) {
        (OrphanMode::Discard, _) if orphan_files => {
            return Err("--orphan1/--orphan2 require --orphan-mode orphan".into());
        }
        (OrphanMode::KeepBoth, _) if orphan_files => {
            return Err("--orphan1/--orphan2 require --orphan-mode orphan".into());
        }
        (OrphanMode::Orphan, true) if args.orphan1.is_none() => {
            return Err("--orphan-mode orphan requires --orphan1".into());
        }
        (OrphanMode::KeepBoth | OrphanMode::Orphan, false) => {
            return Err(
                "--orphan-mode keep-both/orphan only applies to paired-end input (--r2)".into(),
            );
        }
        _ => {}
    }
    let mut outputs: Vec<(&str, &Path)> = vec![("--o1", &args.o1)];
    if let Some(p) = &args.o2 {
        outputs.push(("--o2", p));
    }
    if let Some(p) = &args.orphan1 {
        outputs.push(("--orphan1", p));
    }
    if let Some(p) = &args.orphan2 {
        outputs.push(("--orphan2", p));
    }
    if let Some(p) = &args.report {
        outputs.push(("--report", p));
    }
    io::check_outputs(&outputs, &inputs, args.force)
}

/// Builds the settings, including adapter auto-detection and the poly-G default.
fn resolve_settings(args: &ShortArgs, log: &Logger) -> Result<Settings> {
    let paired = args.r2.is_some();
    let front2 = args.trim_front2.unwrap_or(args.trim_front1);
    let tail2 = args.trim_tail2.unwrap_or(args.trim_tail1);
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
        right: window(
            args.cut_right,
            args.cut_right_window_size,
            args.cut_right_mean_quality,
            args,
        ),
    };
    if args.cut_right && args.cut_tail {
        log.warn("--cut-right overrides --cut-tail; --cut-tail has no effect");
    }

    let first_name = peek_first_name(&args.r1, "R1")?;
    let header = first_name.as_deref().map(parse_header);
    if let Some(h) = &header
        && (h.platform == Platform::OxfordNanopore || h.platform == Platform::PacBio)
    {
        log.warn(&format!(
            "R1 looks like {} long reads; 'short' is meant for Illumina-style reads \
                 ('simpletrim long' is not implemented yet)",
            h.platform.label()
        ));
    }

    let (poly_g, poly_g_reason) = if args.trim_poly_g {
        (
            Some(args.poly_g_min_len as usize),
            "enabled by --trim-poly-g",
        )
    } else if args.disable_trim_poly_g {
        (None, "disabled by --disable-trim-poly-g")
    } else if header.as_ref().is_some_and(|h| h.two_color) {
        (
            Some(args.poly_g_min_len as usize),
            "auto: two-colour instrument",
        )
    } else {
        (None, "auto: not a two-colour instrument")
    };

    let adapter = if args.disable_adapter_trimming {
        None
    } else {
        Some(resolve_adapters(args, paired, log)?)
    };

    let quality_filter = (!args.disable_quality_filter).then_some(QualityFilter {
        qualified_quality: args.qualified_quality,
        unqualified_percent_limit: args.unqualified_percent_limit,
        n_base_limit: args.n_base_limit,
        average_quality: args.average_quality,
    });
    let length_filter = (!args.disable_length_filter).then_some(LengthFilter {
        min_len: args.min_length,
        max_len: args.max_length,
    });
    if args.max_length > 0 && args.max_length < args.min_length && length_filter.is_some() {
        return Err("--max-length is smaller than --min-length".into());
    }

    let orphan2_slot = if args.orphan2.is_some() {
        ORPHAN2
    } else {
        ORPHAN1
    };
    let gz_of = |p: &Option<PathBuf>| p.as_deref().is_some_and(is_gz);
    let gz = [
        is_gz(&args.o1),
        gz_of(&args.o2),
        gz_of(&args.orphan1),
        gz_of(&args.orphan2),
    ];

    Ok(Settings {
        paired,
        front: [args.trim_front1, front2],
        tail: [args.trim_tail1, tail2],
        cut,
        poly_g,
        poly_g_reason,
        dimer_max_len: (adapter.is_some() && !args.disable_dimer_check)
            .then_some(args.dimer_max_len),
        adapter,
        quality_filter,
        length_filter,
        orphan_mode: args.orphan_mode,
        check_names: !args.no_name_check,
        orphan2_slot,
        gz,
        compression_level: args.compression_level,
    })
}

fn resolve_adapters(args: &ShortArgs, paired: bool, log: &Logger) -> Result<AdapterSettings> {
    let mut a = [
        parse_adapter_arg(&args.adapter_r1, "--adapter-r1")?,
        parse_adapter_arg(&args.adapter_r2, "--adapter-r2")?,
    ];
    if !paired && a[1] != AdapterArg::Auto {
        log.warn("--adapter-r2 has no effect for single-end input");
    }
    if let Some((r1, r2)) = args
        .adapter_preset
        .as_deref()
        .and_then(known_adapters::preset)
    {
        for (arg, seq) in a.iter_mut().zip([r1, r2]) {
            if *arg == AdapterArg::Auto {
                *arg = AdapterArg::Seq(seq.as_bytes().to_vec(), Source::Preset);
            }
        }
    }
    // as in fastp: read 2 uses the read 1 adapter unless it is detected separately
    if paired
        && a[1] == AdapterArg::Auto
        && !args.detect_adapter_pe
        && let AdapterArg::Seq(seq, _) = &a[0]
    {
        a[1] = AdapterArg::Seq(seq.clone(), Source::SameAsR1);
    }

    let detect_mate = |mate: usize| -> bool {
        a[mate] == AdapterArg::Auto
            && (if paired {
                args.detect_adapter_pe
            } else {
                mate == 0
            })
    };
    let paths = [Some(args.r1.as_path()), args.r2.as_deref()];
    let tails = [args.trim_tail1, args.trim_tail2.unwrap_or(args.trim_tail1)];
    let mut detected: [Option<Result<Option<Detected>>>; 2] = [None, None];
    std::thread::scope(|s| {
        let handles: Vec<_> = (0..2)
            .filter(|&m| detect_mate(m) && paths[m].is_some())
            .map(|m| {
                let path = paths[m].unwrap();
                log.info(&format!("Detecting adapter for R{} ...", m + 1));
                (
                    m,
                    s.spawn(move || {
                        detect::load_sample(path, &format!("R{}", m + 1))
                            .map(|seqs| detect::detect_adapter(&seqs, tails[m]))
                    }),
                )
            })
            .collect();
        for (m, h) in handles {
            detected[m] = Some(h.join().expect("adapter detection panicked"));
        }
    });

    let mut seq: [Option<Vec<u8>>; 2] = [None, None];
    let mut source = [None, None];
    let mut names = [None, None];
    for m in 0..2 {
        if let Some(result) = detected[m].take() {
            match result? {
                Some(d) => {
                    log.info_with(
                        "adapter",
                        &format!(
                            "R{} adapter detected: {}{}",
                            m + 1,
                            d.sequence,
                            d.name.map(|n| format!(" ({n})")).unwrap_or_default()
                        ),
                        json!({ "read": format!("r{}", m + 1), "sequence": d.sequence, "name": d.name }),
                    );
                    names[m] = d.name;
                    seq[m] = Some(d.sequence.into_bytes());
                    source[m] = Some(Source::Detected);
                }
                None => log.info_with(
                    "adapter",
                    &format!("R{}: no adapter detected", m + 1),
                    json!({ "read": format!("r{}", m + 1), "sequence": null }),
                ),
            }
        } else if let AdapterArg::Seq(s, src) = &a[m]
            && (m == 0 || paired)
        {
            names[m] = known_adapters::describe(s);
            seq[m] = Some(s.clone());
            source[m] = Some(*src);
        }
    }

    let fasta = match &args.adapter_fasta {
        Some(p) => load_adapter_fasta(p, log)?,
        None => Vec::new(),
    };
    Ok(AdapterSettings {
        seq,
        source,
        names,
        fasta,
        matching: MatchParams {
            min_overlap: args.adapter_min_overlap.map(usize::from),
            error_rate: args.adapter_error_rate,
            max_mismatches: args.adapter_max_mismatches.map(usize::from),
        },
        overlap: OverlapParams {
            min_len: args.overlap_min_len as usize,
            max_mismatches: args.overlap_max_mismatches as usize,
            max_error_rate: args.overlap_error_rate,
            allow_gap: args.allow_gap_overlap_trimming,
        },
    })
}

/// Fixed trimming, sliding-window cut and poly-G trimming of one read.
/// Returns the kept region `[start, end)`, or `None` if nothing is left.
fn trim_ends(
    s: &Settings,
    mate: usize,
    rec: &Record,
    st: &mut TrimStats,
) -> Option<(usize, usize)> {
    let len = rec.seq.len();
    let kept = trim_and_cut(rec.seq, rec.qual, s.front[mate], s.tail[mate], &s.cut);
    let fixed = (s.front[mate] + s.tail[mate]).min(len);
    st.fixed.add(fixed);
    st.quality_window
        .add((len - fixed).saturating_sub(kept.map_or(0, |(a, b)| b - a)));
    let (start, mut end) = kept?;
    if let Some(min_len) = s.poly_g
        && let Some(new_len) = poly_g_trim_len(&rec.seq[start..end], min_len)
    {
        st.poly_g.add(end - start - new_len);
        end = start + new_len;
    }
    Some((start, end))
}

/// Trims the configured adapter sequence and FASTA adapters from one read region.
fn trim_adapter_sequences(
    ad: &AdapterSettings,
    mate: usize,
    seq: &[u8],
    region: &mut (usize, usize),
    use_sequence: bool,
) -> bool {
    let mut trimmed = false;
    if use_sequence
        && let Some(adapter) = &ad.seq[mate]
        && let Some(n) = trim_single(&seq[region.0..region.1], adapter, &ad.matching)
    {
        region.1 = region.0 + n;
        trimmed = true;
    }
    if !ad.fasta.is_empty()
        && let Some(n) = trim_multi(&seq[region.0..region.1], &ad.fasta, &ad.matching)
    {
        region.1 = region.0 + n;
        trimmed = true;
    }
    trimmed
}

fn fate_of(s: &Settings, rec: &Record, kept: Option<(usize, usize)>) -> Fate {
    match kept {
        None => Fate::TooShort,
        Some((a, b)) => filter::evaluate(
            &rec.seq[a..b],
            &rec.qual[a..b],
            s.quality_filter.as_ref(),
            s.length_filter.as_ref(),
        ),
    }
}

fn process_single(s: &Settings, rec: &Record, st: &mut RunStats, out: &mut [Vec<u8>; SLOTS]) {
    st.input[0].add(rec.seq, rec.qual);
    let mut kept = trim_ends(s, 0, rec, &mut st.trim[0]);
    let mut dimer = false;
    if let (Some(ad), Some(region)) = (&s.adapter, kept.as_mut()) {
        let before = region.1 - region.0;
        if trim_adapter_sequences(ad, 0, rec.seq, region, true) {
            st.trim[0].adapter.reads += 1;
            st.trim[0].adapter.bases += (before - (region.1 - region.0)) as u64;
            dimer = s
                .dimer_max_len
                .is_some_and(|max| region.1 - region.0 <= max);
        }
    }
    let fate = if dimer {
        Fate::AdapterDimer
    } else {
        fate_of(s, rec, kept)
    };
    st.fates[0].add(fate);
    if fate == Fate::Pass {
        let (a, b) = kept.expect("passing reads are not empty");
        rec.write_trimmed(a, b, &mut out[OUT1]);
        st.output[0].add(&rec.seq[a..b], &rec.qual[a..b]);
    }
}

fn process_pair(
    s: &Settings,
    recs: [&Record; 2],
    st: &mut RunStats,
    out: &mut [Vec<u8>; SLOTS],
    rc_buf: &mut Vec<u8>,
) {
    st.pairs.pairs += 1;
    let mut kept = [None, None];
    for m in 0..2 {
        st.input[m].add(recs[m].seq, recs[m].qual);
        kept[m] = trim_ends(s, m, recs[m], &mut st.trim[m]);
    }

    let mut dimer = false;
    if let Some(ad) = &s.adapter {
        let before = kept.map(|k| k.map_or(0, |(a, b)| b - a));
        let mut trimmed = [false; 2];
        let mut by_overlap = false;
        if let [Some(k1), Some(k2)] = &mut kept {
            let ov = overlap::analyze(
                &recs[0].seq[k1.0..k1.1],
                &recs[1].seq[k2.0..k2.1],
                &ad.overlap,
                rc_buf,
            );
            if let Some((n1, n2)) = overlap::adapter_cut(k1.1 - k1.0, k2.1 - k2.0, &ov, k1.0, k2.0)
            {
                k1.1 = k1.0 + n1;
                k2.1 = k2.0 + n2;
                trimmed = [true, true];
                by_overlap = true;
                st.pairs.adapter_trimmed_by_overlap += 1;
            }
        }
        for m in 0..2 {
            if let Some(region) = kept[m].as_mut() {
                trimmed[m] |= trim_adapter_sequences(ad, m, recs[m].seq, region, !by_overlap);
            }
            if trimmed[m] {
                let after = kept[m].map_or(0, |(a, b)| b - a);
                st.trim[m].adapter.reads += 1;
                st.trim[m].adapter.bases += (before[m] - after) as u64;
            }
        }
        if let (Some(max), [Some(k1), Some(k2)]) = (s.dimer_max_len, kept) {
            dimer = (trimmed[0] || trimmed[1]) && k1.1 - k1.0 <= max && k2.1 - k2.0 <= max;
        }
    }

    let fates = if dimer {
        st.pairs.adapter_dimers += 1;
        [Fate::AdapterDimer; 2]
    } else {
        [fate_of(s, recs[0], kept[0]), fate_of(s, recs[1], kept[1])]
    };
    st.fates[0].add(fates[0]);
    st.fates[1].add(fates[1]);

    let write_trimmed = |m: usize, slot: usize, out: &mut [Vec<u8>; SLOTS]| {
        let (a, b) = kept[m].expect("passing reads are not empty");
        recs[m].write_trimmed(a, b, &mut out[slot]);
        (a, b)
    };
    match (fates[0] == Fate::Pass, fates[1] == Fate::Pass) {
        (true, true) => {
            st.pairs.both_passed += 1;
            st.pairs.pairs_written += 1;
            for (m, rec) in recs.iter().enumerate() {
                let (a, b) = write_trimmed(m, m, out);
                st.output[m].add(&rec.seq[a..b], &rec.qual[a..b]);
            }
        }
        (false, false) => st.pairs.neither_passed += 1,
        (p1, _) => {
            let (pass, fail) = if p1 { (0, 1) } else { (1, 0) };
            if p1 {
                st.pairs.only_r1_passed += 1;
            } else {
                st.pairs.only_r2_passed += 1;
            }
            match s.orphan_mode {
                OrphanMode::Discard => {}
                OrphanMode::KeepBoth => {
                    st.pairs.pairs_written += 1;
                    let (a, b) = write_trimmed(pass, pass, out);
                    st.output[pass].add(&recs[pass].seq[a..b], &recs[pass].qual[a..b]);
                    recs[fail].write(&mut out[fail]);
                    st.output[fail].add(recs[fail].seq, recs[fail].qual);
                    if fail == 0 {
                        st.pairs.rescued_r1_untrimmed += 1;
                    } else {
                        st.pairs.rescued_r2_untrimmed += 1;
                    }
                }
                OrphanMode::Orphan => {
                    let slot = if pass == 0 { ORPHAN1 } else { s.orphan2_slot };
                    let (a, b) = write_trimmed(pass, slot, out);
                    st.orphans
                        .add(&recs[pass].seq[a..b], &recs[pass].qual[a..b]);
                    if pass == 0 {
                        st.pairs.orphans_r1_written += 1;
                    } else {
                        st.pairs.orphans_r2_written += 1;
                    }
                }
            }
        }
    }
}

struct Job {
    index: usize,
    r1: Chunk,
    r2: Option<Chunk>,
}

struct Done {
    index: usize,
    records: u64,
    outputs: [Vec<u8>; SLOTS],
    stats: RunStats,
}

fn process_job(s: &Settings, job: Job) -> Result<Done> {
    let mut st = RunStats::default();
    let mut out: [Vec<u8>; SLOTS] = Default::default();
    match &job.r2 {
        None => {
            for rec in job.r1.iter() {
                process_single(s, &rec, &mut st, &mut out);
            }
        }
        Some(r2) => {
            let mut rc_buf = Vec::new();
            for i in 0..job.r1.len() {
                let (a, b) = (job.r1.get(i), r2.get(i));
                if s.check_names && !mate_names_match(&a, &b) {
                    return Err(format!(
                        "read names of pair {} do not match: '{}' (R1) vs '{}' (R2); \
                         the files may be out of sync (use --no-name-check to skip this check)",
                        job.r1.first_record + i as u64 + 1,
                        String::from_utf8_lossy(a.id()),
                        String::from_utf8_lossy(b.id()),
                    ));
                }
                process_pair(s, [&a, &b], &mut st, &mut out, &mut rc_buf);
            }
        }
    }
    for (slot, data) in out.iter_mut().enumerate() {
        if s.gz[slot] && !data.is_empty() {
            *data = gzip_member(data, s.compression_level);
        }
    }
    Ok(Done {
        index: job.index,
        records: job.r1.len() as u64,
        outputs: out,
        stats: st,
    })
}

/// An output file being written.
struct Sink {
    file: File,
    path: PathBuf,
    gz: bool,
    written: u64,
}

impl Sink {
    fn write(&mut self, data: &[u8]) -> Result<()> {
        self.file
            .write_all(data)
            .map_err(|e| format!("cannot write {}: {e}", self.path.display()))?;
        self.written += data.len() as u64;
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        if self.gz && self.written == 0 {
            // an empty file is not valid gzip; write an empty member instead
            self.write(&gzip_member(b"", 1))?;
        }
        self.file
            .flush()
            .map_err(|e| format!("cannot write {}: {e}", self.path.display()))
    }
}

/// Periodic progress reports from the writer thread.
struct Progress<'a> {
    log: &'a Logger,
    interval: Option<Duration>,
    start: Instant,
    last: Instant,
    unit: &'static str,
    bytes_read: Arc<AtomicU64>,
    file_size: u64,
}

impl Progress<'_> {
    fn update(&mut self, records: u64) {
        let Some(interval) = self.interval else {
            return;
        };
        if self.last.elapsed() < interval {
            return;
        }
        self.last = Instant::now();
        let secs = self.start.elapsed().as_secs_f64().max(1e-9);
        let rate = records as f64 / secs;
        let fraction = (self.file_size > 0).then(|| {
            (self.bytes_read.load(Ordering::Relaxed) as f64 / self.file_size as f64).min(1.0)
        });
        let pct = fraction
            .map(|f| format!(" ({:.1}% of input)", f * 100.0))
            .unwrap_or_default();
        self.log.progress(
            &format!(
                "processed {} {}{pct}, {} {}/s",
                group(records),
                self.unit,
                group(rate as u64),
                self.unit
            ),
            json!({
                "records": records,
                "unit": self.unit,
                "input_fraction": fraction.map(|f| (f * 10000.0).round() / 10000.0),
                "records_per_s": rate.round(),
            }),
        );
    }
}

fn write_loop(
    rx: Receiver<Result<Done>>,
    tokens: SyncSender<()>,
    mut sinks: Vec<Option<Sink>>,
    mut progress: Progress,
) -> Result<(RunStats, u64)> {
    let mut pending: BTreeMap<usize, Done> = BTreeMap::new();
    let mut next = 0;
    let mut stats = RunStats::default();
    let mut records = 0u64;
    for msg in rx {
        let done = msg?;
        pending.insert(done.index, done);
        while let Some(done) = pending.remove(&next) {
            for (slot, data) in done.outputs.iter().enumerate() {
                if !data.is_empty() {
                    sinks[slot]
                        .as_mut()
                        .expect("output slot without file")
                        .write(data)?;
                }
            }
            stats.merge(&done.stats);
            records += done.records;
            next += 1;
            let _ = tokens.send(());
            progress.update(records);
        }
    }
    for sink in sinks.into_iter().flatten() {
        sink.finish()?;
    }
    Ok((stats, records))
}

fn receive(rx: Option<&Receiver<Result<Chunk>>>) -> Result<Option<Chunk>> {
    match rx.map(|rx| rx.recv()) {
        Some(Ok(chunk)) => chunk.map(Some),
        _ => Ok(None),
    }
}

/// Runs the reader → workers → writer pipeline.
fn run_pipeline(
    s: &Settings,
    args: &ShortArgs,
    threads: usize,
    sinks: Vec<Option<Sink>>,
    log: &Logger,
) -> Result<(RunStats, u64)> {
    let in1 = io::open_fastq(&args.r1, "R1")?;
    let in2 = args
        .r2
        .as_deref()
        .map(|p| io::open_fastq(p, "R2"))
        .transpose()?;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .thread_name(|i| format!("simpletrim-worker-{i}"))
        .build()
        .map_err(|e| format!("cannot start worker threads: {e}"))?;
    let max_in_flight = threads * 2 + 2;
    let progress = Progress {
        log,
        // zero, negative, NaN or out-of-range values disable progress reports
        interval: Duration::try_from_secs_f64(args.progress_interval)
            .ok()
            .filter(|d| !d.is_zero()),
        start: Instant::now(),
        last: Instant::now(),
        unit: if s.paired { "pairs" } else { "reads" },
        bytes_read: in1.bytes_read.clone(),
        file_size: in1.file_size,
    };

    std::thread::scope(|scope| {
        let (tx1, rx1) = sync_channel(max_in_flight);
        let r1_thread =
            scope.spawn(move || feed_chunks(in1.reader, tx1, args.max_reads, CHUNK_RECORDS));
        let (rx2, r2_thread) = match in2 {
            Some(in2) => {
                let (tx2, rx2) = sync_channel(max_in_flight);
                let t = scope
                    .spawn(move || feed_chunks(in2.reader, tx2, args.max_reads, CHUNK_RECORDS));
                (Some(rx2), Some(t))
            }
            None => (None, None),
        };
        let (done_tx, done_rx) = channel::<Result<Done>>();
        let (token_tx, token_rx) = sync_channel::<()>(max_in_flight);
        for _ in 0..max_in_flight {
            token_tx.send(()).expect("token channel has capacity");
        }
        let writer = scope.spawn(move || write_loop(done_rx, token_tx, sinks, progress));

        let mut dispatch_error = None;
        pool.in_place_scope(|ps| {
            let mut index = 0;
            loop {
                if token_rx.recv().is_err() {
                    break; // the writer stopped (with an error)
                }
                let next = receive(Some(&rx1)).and_then(|c1| Ok((c1, receive(rx2.as_ref())?)));
                let job = match next {
                    Err(e) => Err(e),
                    Ok((None, None)) => break,
                    Ok((Some(r1), None)) if !s.paired => Ok(Job { index, r1, r2: None }),
                    Ok((Some(r1), Some(r2))) if r1.len() == r2.len() => Ok(Job { index, r1, r2: Some(r2) }),
                    Ok((c1, c2)) => {
                        let reads = |c: &Option<Chunk>, other: &Option<Chunk>| match (c, other) {
                            (Some(c), _) => c.first_record + c.len() as u64,
                            (None, Some(o)) => o.first_record,
                            (None, None) => 0,
                        };
                        let (n1, n2) = (reads(&c1, &c2), reads(&c2, &c1));
                        let (shorter, n) = if n1 < n2 { ("R1", n1) } else { ("R2", n2) };
                        Err(format!(
                            "R1 and R2 have different numbers of reads: {shorter} ended after {} reads",
                            group(n)
                        ))
                    }
                };
                match job {
                    Ok(job) => {
                        let tx = done_tx.clone();
                        ps.spawn(move |_| {
                            let _ = tx.send(process_job(s, job));
                        });
                        index += 1;
                    }
                    Err(e) => {
                        dispatch_error = Some(e);
                        break;
                    }
                }
            }
        });
        drop(done_tx);
        drop(rx1);
        drop(rx2);
        let written = writer.join().expect("writer thread panicked");
        r1_thread.join().expect("reader thread panicked");
        if let Some(t) = r2_thread {
            t.join().expect("reader thread panicked");
        }
        match dispatch_error {
            Some(e) => Err(e),
            None => written,
        }
    })
}

fn describe_settings(s: &Settings) -> (Vec<String>, Value) {
    let mut lines = Vec::new();
    let fixed = s.front.iter().chain(&s.tail).any(|&n| n > 0);
    lines.push(if fixed {
        if s.paired {
            format!(
                "fixed trimming: R1 front {} tail {}, R2 front {} tail {}",
                s.front[0], s.tail[0], s.front[1], s.tail[1]
            )
        } else {
            format!("fixed trimming: front {} tail {}", s.front[0], s.tail[0])
        }
    } else {
        "fixed trimming: off".into()
    });
    let win = |name: &str, w: Option<Window>| {
        w.map(|w| format!("{name} (window {}, Q{})", w.size, w.quality))
    };
    let cuts: Vec<String> = [
        win("front", s.cut.front),
        win(
            "tail",
            if s.cut.right.is_some() {
                None
            } else {
                s.cut.tail
            },
        ),
        win("right", s.cut.right),
    ]
    .into_iter()
    .flatten()
    .collect();
    lines.push(format!(
        "quality window: {}",
        if cuts.is_empty() {
            "off".into()
        } else {
            cuts.join(", ")
        }
    ));
    lines.push(match s.poly_g {
        Some(n) => format!("poly-G trimming: on, min length {n} ({})", s.poly_g_reason),
        None => format!("poly-G trimming: off ({})", s.poly_g_reason),
    });
    let mut adapter_json = Value::Null;
    match &s.adapter {
        None => lines.push("adapter trimming: off".into()),
        Some(ad) => {
            let mates = if s.paired { 2 } else { 1 };
            let mut parts = Vec::new();
            if s.paired {
                parts.push("read overlap".to_string());
            }
            for m in 0..mates {
                if let Some(seq) = &ad.seq[m] {
                    parts.push(format!(
                        "R{} {} ({}{})",
                        m + 1,
                        String::from_utf8_lossy(seq),
                        ad.source[m].map_or("", Source::name),
                        ad.names[m].map(|n| format!(": {n}")).unwrap_or_default()
                    ));
                }
            }
            if !ad.fasta.is_empty() {
                parts.push(format!("{} FASTA sequences", ad.fasta.len()));
            }
            if parts.is_empty() {
                parts.push("no adapter sequence".into());
            }
            lines.push(format!("adapter trimming: {}", parts.join(" + ")));
            let mate_json = |m: usize| {
                json!({
                    "sequence": ad.seq[m].as_ref().map(|s| String::from_utf8_lossy(s).into_owned()),
                    "source": ad.source[m].map(Source::name),
                    "name": ad.names[m],
                })
            };
            adapter_json = json!({
                "r1": mate_json(0),
                "r2": if s.paired { mate_json(1) } else { Value::Null },
                "fasta_sequences": ad.fasta.len(),
                "overlap": if s.paired { json!({
                    "min_len": ad.overlap.min_len,
                    "max_mismatches": ad.overlap.max_mismatches,
                    "max_error_rate": ad.overlap.max_error_rate,
                    "allow_gap": ad.overlap.allow_gap,
                }) } else { Value::Null },
                "matching": {
                    "min_overlap": ad.matching.min_overlap,
                    "error_rate": ad.matching.error_rate,
                    "max_mismatches": ad.matching.max_mismatches,
                },
            });
        }
    }
    lines.push(match s.dimer_max_len {
        Some(n) => format!("adapter dimer check: on (<= {n} bp after adapter trimming)"),
        None => "adapter dimer check: off".into(),
    });
    lines.push(match &s.quality_filter {
        Some(q) => format!(
            "quality filter: <= {}% bases below Q{}, <= {} N{}",
            q.unqualified_percent_limit,
            q.qualified_quality,
            q.n_base_limit,
            if q.average_quality > 0 {
                format!(", mean >= Q{}", q.average_quality)
            } else {
                String::new()
            }
        ),
        None => "quality filter: off".into(),
    });
    lines.push(match &s.length_filter {
        Some(l) if l.max_len > 0 => format!("length filter: {}-{} bp", l.min_len, l.max_len),
        Some(l) => format!("length filter: >= {} bp", l.min_len),
        None => "length filter: off".into(),
    });
    if s.paired {
        lines.push(format!(
            "pairs with one failing read: {}{}",
            s.orphan_mode.name(),
            if s.check_names {
                ""
            } else {
                " (name check off)"
            }
        ));
    }
    let json = json!({
        "trim_front": if s.paired { json!(s.front) } else { json!([s.front[0]]) },
        "trim_tail": if s.paired { json!(s.tail) } else { json!([s.tail[0]]) },
        "cut_front": s.cut.front.map(|w| json!({"window": w.size, "quality": w.quality})),
        "cut_tail": s.cut.tail.filter(|_| s.cut.right.is_none()).map(|w| json!({"window": w.size, "quality": w.quality})),
        "cut_right": s.cut.right.map(|w| json!({"window": w.size, "quality": w.quality})),
        "poly_g": { "enabled": s.poly_g.is_some(), "min_len": s.poly_g, "reason": s.poly_g_reason },
        "adapter": adapter_json,
        "dimer_max_len": s.dimer_max_len,
        "quality_filter": s.quality_filter.map(|q| json!({
            "qualified_quality": q.qualified_quality,
            "unqualified_percent_limit": q.unqualified_percent_limit,
            "n_base_limit": q.n_base_limit,
            "average_quality": q.average_quality,
        })),
        "length_filter": s.length_filter.map(|l| json!({"min_len": l.min_len, "max_len": l.max_len})),
        "orphan_mode": if s.paired { Some(s.orphan_mode.name()) } else { None },
        "check_names": s.paired && s.check_names,
        "compression_level": s.compression_level,
    });
    (lines, json)
}

fn summary(s: &Settings, args: &ShortArgs, st: &RunStats, elapsed: f64) -> (String, Value) {
    let mates = if s.paired { 2 } else { 1 };
    let unit = if s.paired { "read pairs" } else { "reads" };
    let total = st.input[0].reads;
    let written = if s.paired {
        st.pairs.pairs_written
    } else {
        st.output[0].reads
    };
    let mut t = String::new();
    let line = |t: &mut String, l: String| {
        t.push_str(&l);
        t.push('\n');
    };
    line(&mut t, format!("Done in {elapsed:.1} s."));
    line(&mut t, format!("Input:  {} {unit}", group(total)));
    let base_line = |b: &crate::stats::BaseStats| {
        format!(
            "{} bases, mean length {:.1}, Q30 {:.2}%, GC {:.2}%",
            group(b.bases),
            b.mean_length(),
            b.q30_rate() * 100.0,
            b.gc_content() * 100.0
        )
    };
    for m in 0..mates {
        line(&mut t, format!("  R{}: {}", m + 1, base_line(&st.input[m])));
    }
    line(
        &mut t,
        format!(
            "Output: {} {unit} ({})",
            group(written),
            percent(written, total)
        ),
    );
    for m in 0..mates {
        line(
            &mut t,
            format!("  R{}: {}", m + 1, base_line(&st.output[m])),
        );
    }
    if st.orphans.reads > 0 || s.orphan_mode == OrphanMode::Orphan {
        line(
            &mut t,
            format!(
                "  orphans: {} (R1 {}, R2 {})",
                group(st.orphans.reads),
                group(st.pairs.orphans_r1_written),
                group(st.pairs.orphans_r2_written)
            ),
        );
    }
    let pair_cell = |a: u64, b: u64| {
        if s.paired {
            format!("{:>14} {:>14}", group(a), group(b))
        } else {
            format!("{:>14}", group(a))
        }
    };
    line(
        &mut t,
        if s.paired {
            format!("Read fates:        {:>14} {:>14}", "R1", "R2")
        } else {
            "Read fates:".into()
        },
    );
    let f = &st.fates;
    for (name, a, b) in [
        ("passed", f[0].passed, f[1].passed),
        ("too short", f[0].too_short, f[1].too_short),
        ("too long", f[0].too_long, f[1].too_long),
        ("low quality", f[0].low_quality, f[1].low_quality),
        ("too many N", f[0].too_many_n, f[1].too_many_n),
        ("adapter dimer", f[0].adapter_dimer, f[1].adapter_dimer),
    ] {
        line(&mut t, format!("  {name:<16} {}", pair_cell(a, b)));
    }
    line(
        &mut t,
        if s.paired {
            "Trimmed (reads / bases):         R1                            R2".into()
        } else {
            "Trimmed (reads / bases):".into()
        },
    );
    let tr = &st.trim;
    let cell =
        |x: &crate::stats::Trimmed| format!("{:>12} / {:<15}", group(x.reads), group(x.bases));
    for (name, a, b) in [
        ("fixed", &tr[0].fixed, &tr[1].fixed),
        (
            "quality window",
            &tr[0].quality_window,
            &tr[1].quality_window,
        ),
        ("poly-G", &tr[0].poly_g, &tr[1].poly_g),
        ("adapter", &tr[0].adapter, &tr[1].adapter),
    ] {
        if s.paired {
            line(&mut t, format!("  {name:<16} {} {}", cell(a), cell(b)));
        } else {
            line(&mut t, format!("  {name:<16} {}", cell(a)));
        }
    }
    if s.paired {
        let p = &st.pairs;
        line(
            &mut t,
            format!(
                "Pairs: both passed {}, only R1 {}, only R2 {}, neither {} (adapter dimers {}); \
                 adapters found by overlap in {} pairs",
                group(p.both_passed),
                group(p.only_r1_passed),
                group(p.only_r2_passed),
                group(p.neither_passed),
                group(p.adapter_dimers),
                group(p.adapter_trimmed_by_overlap)
            ),
        );
        if s.orphan_mode == OrphanMode::KeepBoth {
            line(
                &mut t,
                format!(
                    "  kept untrimmed because the mate passed: R1 {}, R2 {}",
                    group(p.rescued_r1_untrimmed),
                    group(p.rescued_r2_untrimmed)
                ),
            );
        }
    }

    let path = |p: &Option<PathBuf>| p.as_ref().map(|p| p.display().to_string());
    let per_mate = |f: &dyn Fn(usize) -> Value| {
        if s.paired {
            json!({"r1": f(0), "r2": f(1)})
        } else {
            json!({"r1": f(0)})
        }
    };
    let json = json!({
        "mode": if s.paired { "paired" } else { "single" },
        "elapsed_s": crate::stats::round(elapsed, 3),
        "inputs": { "r1": args.r1.display().to_string(), "r2": path(&args.r2) },
        "outputs": {
            "o1": args.o1.display().to_string(),
            "o2": path(&args.o2),
            "orphan1": path(&args.orphan1),
            "orphan2": path(&args.orphan2),
        },
        "records_in": total,
        "records_written": written,
        "before": per_mate(&|m| st.input[m].to_json()),
        "after": per_mate(&|m| st.output[m].to_json()),
        "orphans": if s.paired { st.orphans.to_json() } else { Value::Null },
        "trimming": per_mate(&|m| serde_json::to_value(st.trim[m]).expect("serializable")),
        "filtering": per_mate(&|m| serde_json::to_value(st.fates[m]).expect("serializable")),
        "pairs": if s.paired { serde_json::to_value(st.pairs).expect("serializable") } else { Value::Null },
    });
    (t, json)
}

/// Runs the `short` subcommand.
pub fn run(args: &ShortArgs, threads: usize, log: &Logger) -> Result<()> {
    validate_paths(args)?;
    let settings = resolve_settings(args, log)?;
    let (lines, settings_json) = describe_settings(&settings);
    let mode = if settings.paired {
        "paired-end"
    } else {
        "single-end"
    };
    let mut header = format!(
        "simpletrim {} short | {mode} | {threads} worker threads\n",
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
            "mode": if settings.paired { "paired" } else { "single" },
            "threads": threads,
            "settings": settings_json,
        }),
    );

    let mut files = OutputFiles::new();
    let mut sinks: Vec<Option<Sink>> = (0..SLOTS).map(|_| None).collect();
    let outputs = [
        (OUT1, Some(&args.o1)),
        (OUT2, args.o2.as_ref()),
        (ORPHAN1, args.orphan1.as_ref()),
        (ORPHAN2, args.orphan2.as_ref()),
    ];
    for (slot, path) in outputs {
        if let Some(path) = path {
            sinks[slot] = Some(Sink {
                file: files.create(path, args.force)?,
                path: path.clone(),
                gz: settings.gz[slot],
                written: 0,
            });
        }
    }

    let (stats, _) = run_pipeline(&settings, args, threads, sinks, log)?;
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
