//! `analyze` for long reads (Oxford Nanopore, PacBio): read lengths and qualities, quality
//! along the reads, adapters at the read ends and inside reads (with the detection `long`
//! uses), poly-X tails and low-quality stretches; suggests a `simpletrim long` command.
//!
//! Quality means average error probabilities, as in `long`.

use std::sync::atomic::Ordering;
use std::sync::mpsc::sync_channel;

use rayon::prelude::*;
use serde_json::{Value, json};

use super::{Suggestion, emit, file_json, header_json, shell_quote, trimmed_name};
use crate::Result;
use crate::cli::AnalyzeArgs;
use crate::fastq::{Chunk, feed_chunks};
use crate::io;
use crate::log::Logger;
use crate::long_adapter::{LongAdapters, Params};
use crate::long_detect::{Detection, Detector};
use crate::platform::HeaderInfo;
use crate::polyx::poly_x_trim_len;
use crate::quality::{Averaging, Scale, error_sum, low_quality_regions, phred_of_error};
use crate::stats::{ReadSetStats, group, round};

const CHUNK_RECORDS: usize = 8192;
const CHUNK_BYTES: usize = 8 << 20;
/// Quality along the reads is reported in bins of this size, from each read end.
const PROFILE_BIN: usize = 25;
const PROFILE_BINS: usize = 10;
/// Reads (and bases) kept for checking where the detected adapters occur.
const ADAPTER_SAMPLE_READS: usize = 5000;
const ADAPTER_SAMPLE_BASES: usize = 50_000_000;
/// Low-quality stretches closer than this to a read end do not count as internal.
const INTERNAL_MARGIN: usize = 150;
/// `long` defaults used for the checks.
const POLY_X_MIN_LEN: usize = 10;
const LOW_QUALITY_WINDOW: usize = 100;
const LOW_QUALITY_THRESHOLD: u8 = 7;
const ADAPTER_PARAMS: Params = Params {
    rate: 0.25,
    min_overlap: 16,
    window: 200,
    extension: 10,
    internal_rate: 0.15,
};
/// A suggestion is made if a feature is seen in at least this fraction of reads.
const SUGGEST_RATE: f64 = 0.05;
/// Yield is reported for reads of at least these lengths.
const YIELD_LENGTHS: [usize; 4] = [1000, 5000, 10_000, 20_000];

#[derive(Clone, Default)]
struct Partial {
    reads: ReadSetStats,
    /// Error-probability sums and base counts per bin, from the read start and from the end.
    from_start: [(f64, u64); PROFILE_BINS],
    from_end: [(f64, u64); PROFILE_BINS],
    /// Reads with a poly-X tail, per tail base (A, C, G, T).
    poly_x: [u64; 4],
    /// Reads with a low-quality stretch away from the read ends, and its bases.
    internal_low_quality: u64,
    internal_low_quality_bases: u64,
    u_reads: u64,
    t_reads: u64,
}

impl Partial {
    fn add(&mut self, seq: &[u8], qual: &[u8]) {
        self.reads.add(seq, qual);
        let len = seq.len();
        for bin in 0..PROFILE_BINS {
            let (a, b) = (bin * PROFILE_BIN, ((bin + 1) * PROFILE_BIN).min(len));
            if a >= b {
                break;
            }
            self.from_start[bin].0 += error_sum(&qual[a..b]);
            self.from_start[bin].1 += (b - a) as u64;
            self.from_end[bin].0 += error_sum(&qual[len - b..len - a]);
            self.from_end[bin].1 += (b - a) as u64;
        }
        if let Some(new_len) = poly_x_trim_len(seq, POLY_X_MIN_LEN) {
            let base = match seq[new_len] {
                b'A' => 0,
                b'C' => 1,
                b'G' => 2,
                _ => 3,
            };
            self.poly_x[base] += 1;
        }
        let errp = Scale::get(Averaging::ErrorProbability);
        let internal: usize =
            low_quality_regions(qual, LOW_QUALITY_WINDOW, LOW_QUALITY_THRESHOLD, errp)
                .into_iter()
                .filter(|&(a, b)| a >= INTERNAL_MARGIN && b + INTERNAL_MARGIN <= len)
                .map(|(a, b)| b - a)
                .sum();
        if internal > 0 {
            self.internal_low_quality += 1;
            self.internal_low_quality_bases += internal as u64;
        }
        self.u_reads += seq.contains(&b'U') as u64;
        self.t_reads += seq.contains(&b'T') as u64;
    }

    fn merge(&mut self, o: &Partial) {
        self.reads.merge(&o.reads);
        for i in 0..PROFILE_BINS {
            self.from_start[i].0 += o.from_start[i].0;
            self.from_start[i].1 += o.from_start[i].1;
            self.from_end[i].0 += o.from_end[i].0;
            self.from_end[i].1 += o.from_end[i].1;
        }
        for i in 0..4 {
            self.poly_x[i] += o.poly_x[i];
        }
        self.internal_low_quality += o.internal_low_quality;
        self.internal_low_quality_bases += o.internal_low_quality_bases;
        self.u_reads += o.u_reads;
        self.t_reads += o.t_reads;
    }
}

fn analyze_chunk(c: &Chunk) -> Partial {
    let mut p = Partial::default();
    for rec in c.iter() {
        p.add(rec.seq, rec.qual);
    }
    p
}

struct Collected {
    stats: Partial,
    detection: Detection,
    sample: Vec<Vec<u8>>,
    whole_file: bool,
    estimated_total: Option<u64>,
}

fn collect(args: &AnalyzeArgs, threads: usize, log: &Logger) -> Result<Collected> {
    let input = io::open_fastq(&args.r1, "input")?;
    let bytes_read = input.bytes_read.clone();
    let file_size = input.file_size;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| format!("cannot start worker threads: {e}"))?;
    let batch_size = threads * 2;
    let max_bases = if args.max_bases == 0 {
        u64::MAX
    } else {
        args.max_bases
    };

    std::thread::scope(|scope| {
        let (tx, rx) = sync_channel(batch_size);
        scope.spawn(move || {
            feed_chunks(input.reader, tx, args.max_reads, CHUNK_RECORDS, CHUNK_BYTES)
        });
        let mut stats = Partial::default();
        let mut detector = Detector::new(1);
        let mut sample = Vec::new();
        let mut sample_bases = 0usize;
        let mut next_report = 100_000u64;
        let mut ended = false;
        while !ended {
            let mut batch = Vec::with_capacity(batch_size);
            while batch.len() < batch_size {
                match rx.recv() {
                    Ok(chunk) => batch.push(chunk?),
                    Err(_) => {
                        ended = true;
                        break;
                    }
                }
            }
            for chunk in &batch {
                for rec in chunk.iter() {
                    if detector.wants_more() {
                        detector.add(rec.seq);
                    }
                    if sample.len() < ADAPTER_SAMPLE_READS && sample_bases < ADAPTER_SAMPLE_BASES {
                        sample_bases += rec.seq.len();
                        sample.push(rec.seq.to_vec());
                    }
                }
            }
            let partials: Vec<Partial> =
                pool.install(|| batch.par_iter().map(analyze_chunk).collect());
            for p in &partials {
                stats.merge(p);
            }
            if stats.reads.base.reads >= next_report {
                log.progress(
                    &format!("read {} reads", group(stats.reads.base.reads)),
                    json!({ "records": stats.reads.base.reads }),
                );
                next_report += 100_000;
            }
            if stats.reads.base.bases >= max_bases {
                break;
            }
        }
        drop(rx);
        let reads = stats.reads.base.reads;
        let consumed = bytes_read.load(Ordering::Relaxed);
        let whole_file = ended && (args.max_reads == 0 || reads < args.max_reads)
            || (file_size > 0 && consumed >= file_size);
        let estimated_total = (!whole_file && consumed > 0)
            .then(|| (reads as f64 * file_size as f64 / consumed as f64).round() as u64);
        Ok(Collected {
            stats,
            detection: detector.finish(),
            sample,
            whole_file,
            estimated_total,
        })
    })
}

/// Fields of nanopore read headers: Dorado (`RG:Z:`, `PU:Z:` ... tags) and MinKNOW
/// (`runid=`, `flow_cell_id=` ...).
fn run_info(name: &[u8]) -> Value {
    let text = String::from_utf8_lossy(name);
    let mut info = serde_json::Map::new();
    let mut set = |key: &str, value: &str| {
        if !value.is_empty() && !info.contains_key(key) {
            info.insert(key.into(), json!(value));
        }
    };
    for token in text.split_ascii_whitespace().skip(1) {
        if let Some((tag, value)) = token
            .split_once(':')
            .and_then(|(t, rest)| rest.split_once(':').map(|(_, v)| (t, v)))
            .filter(|(t, _)| t.len() == 2)
        {
            match tag {
                "RG" => {
                    // <run id>_<basecall model>_<barcode>
                    if let Some((run, rest)) = value.split_once('_') {
                        set("run_id", run);
                        let model_end = rest
                            .find('@')
                            .map(|at| rest[at..].find('_').map_or(rest.len(), |e| at + e))
                            .unwrap_or(rest.len());
                        set("basecall_model", &rest[..model_end]);
                    }
                }
                "PU" => set("flow_cell", value),
                "al" => set("barcode", value),
                "DT" => set("start_time", value),
                _ => {}
            }
        } else if let Some((key, value)) = token.split_once('=') {
            match key {
                "runid" => set("run_id", value),
                "flow_cell_id" => set("flow_cell", value),
                "basecall_model_version_id" | "model_version_id" => set("basecall_model", value),
                "barcode" => set("barcode", value),
                "start_time" => set("start_time", value),
                "sample_id" => set("sample", value),
                _ => {}
            }
        }
    }
    Value::Object(info)
}

/// Fraction of sampled reads with the adapters at the start, at the end and inside.
fn adapter_prevalence(
    sample: &[Vec<u8>],
    adapters: &LongAdapters,
    pool: &rayon::ThreadPool,
) -> [f64; 3] {
    if sample.is_empty() {
        return [0.0; 3];
    }
    let counts = pool.install(|| {
        sample
            .par_iter()
            .map_init(Vec::new, |buf, seq| {
                let a = adapters.trim_start(seq);
                let b = seq.len() - adapters.trim_end(&seq[a..], buf);
                let inside = b > a && !adapters.internal_regions(&seq[a..b]).is_empty();
                [(a > 0) as u64, (b < seq.len()) as u64, inside as u64]
            })
            .reduce(|| [0; 3], |x, y| [x[0] + y[0], x[1] + y[1], x[2] + y[2]])
    });
    counts.map(|n| n as f64 / sample.len() as f64)
}

fn q_or_null(err: f64, n: u64) -> Value {
    if n == 0 {
        Value::Null
    } else {
        json!(round(phred_of_error(err / n as f64), 1))
    }
}

/// Runs `analyze` for long-read input.
pub fn run(
    args: &AnalyzeArgs,
    header: &HeaderInfo,
    first_name: &[u8],
    threads: usize,
    log: &Logger,
) -> Result<()> {
    if args.r2.is_some() {
        return Err("--r2 is not supported for long reads (they are single-end)".into());
    }
    log.info(&format!(
        "Reading {} ...",
        match (args.max_reads, args.max_bases) {
            (0, 0) => "all reads".to_string(),
            (0, b) => format!("up to {} bases", group(b)),
            (r, 0) => format!("up to {} reads", group(r)),
            (r, b) => format!("up to {} reads or {} bases", group(r), group(b)),
        }
    ));
    let c = collect(args, threads, log)?;
    let st = &c.stats;
    let rs = &st.reads;
    let total = rs.base.reads;
    if total == 0 {
        return Err(format!("{} contains no reads", args.r1.display()));
    }
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| format!("cannot start worker threads: {e}"))?;
    let d = &c.detection;
    let to_t = |s: &String| {
        s.bytes()
            .map(|b| if b == b'U' { b'T' } else { b })
            .collect::<Vec<u8>>()
    };
    let start = d.start.adapter.as_ref().map(to_t);
    let end = d.end.adapter.as_ref().map(to_t);
    let any_adapter = start.is_some() || end.is_some();
    let adapters = LongAdapters::new(start.as_deref(), end.as_deref(), &[], ADAPTER_PARAMS);
    let [at_start, at_end, inside] = if any_adapter {
        adapter_prevalence(&c.sample, &adapters, &pool)
    } else {
        [0.0; 3]
    };

    let mut warnings: Vec<String> = Vec::new();
    let mut notes: Vec<String> = Vec::new();
    let mut suggestions: Vec<Suggestion> = Vec::new();
    let rate = |n: u64| n as f64 / total as f64;

    let rna = st.u_reads > 0 && st.t_reads == 0;
    if st.u_reads > 0 && st.t_reads > 0 {
        warnings.push("the reads contain both U and T bases; 'simpletrim long' stops with an error on such input".into());
    } else if rna {
        notes.push("the reads contain U instead of T (direct RNA)".into());
    }

    // adapters
    if any_adapter {
        let flag = |s: &Option<Vec<u8>>| {
            s.as_ref().map_or("none".to_string(), |s| {
                String::from_utf8_lossy(s).into_owned()
            })
        };
        suggestions.push(Suggestion {
            flags: vec![
                "--start-adapter".into(),
                flag(&start),
                "--end-adapter".into(),
                flag(&end),
            ],
            reason: "detected adapters (same as the default auto-detection, but skips the detection pass)".into(),
        });
    } else if d.reads >= crate::long_detect::MIN_READS {
        notes.push("no adapters detected at the read starts or ends; the reads may already be trimmed (e.g. by the basecaller)".into());
    } else {
        notes.push(format!(
            "too few reads for adapter detection (need {})",
            crate::long_detect::MIN_READS
        ));
    }
    for (which, flag, e) in [
        ("start", "--start-adapter", &d.start),
        ("end", "--end-adapter", &d.end),
    ] {
        if let Some(s) = &e.too_short {
            notes.push(format!(
                "frequent sequence {s} at the read {which}s is too short to be used as an adapter automatically; pass {flag} {s} if it is one"
            ));
        }
    }
    if inside >= 0.01 {
        notes.push(format!(
            "{:.2}% of the sampled reads contain an adapter inside (chimeric reads); 'long' splits them at the adapter by default",
            inside * 100.0
        ));
    }

    // quality and length
    let poly_x: u64 = st.poly_x.iter().sum();
    if rate(poly_x) >= SUGGEST_RATE {
        suggestions.push(Suggestion {
            flags: vec!["--trim-poly-x".into()],
            reason: format!(
                "{:.2}% of reads end in a homopolymer tail (>= {POLY_X_MIN_LEN} bases)",
                rate(poly_x) * 100.0
            ),
        });
    }
    let first_bin = phred_of_error(st.from_start[0].0 / st.from_start[0].1.max(1) as f64);
    let last_bin = phred_of_error(st.from_end[0].0 / st.from_end[0].1.max(1) as f64);
    if start.is_none() && first_bin < LOW_QUALITY_THRESHOLD as f64 {
        suggestions.push(Suggestion {
            flags: vec!["--cut-front".into()],
            reason: format!("the first {PROFILE_BIN} bases have a mean quality of Q{first_bin:.1} and no start adapter was found"),
        });
    }
    if end.is_none() && last_bin < LOW_QUALITY_THRESHOLD as f64 {
        suggestions.push(Suggestion {
            flags: vec!["--cut-tail".into()],
            reason: format!("the last {PROFILE_BIN} bases have a mean quality of Q{last_bin:.1} and no end adapter was found"),
        });
    }
    if rate(st.internal_low_quality) >= SUGGEST_RATE {
        suggestions.push(Suggestion {
            flags: vec!["--low-quality-mode".into(), "split".into()],
            reason: format!(
                "{:.2}% of reads have a low-quality stretch inside (window {LOW_QUALITY_WINDOW}, mean below Q{LOW_QUALITY_THRESHOLD})",
                rate(st.internal_low_quality) * 100.0
            ),
        });
    }
    let below_q7 = total - rs.reads_with_quality_at_least(7);
    if rate(below_q7) >= SUGGEST_RATE {
        suggestions.push(Suggestion {
            flags: vec!["--average-quality".into(), "7".into()],
            reason: format!(
                "{:.2}% of reads have a mean quality below Q7",
                rate(below_q7) * 100.0
            ),
        });
    }
    if !c.whole_file {
        notes.push(
            "only the first part of the file was inspected (see --max-reads and --max-bases)"
                .into(),
        );
    }

    let mut cmd = vec![
        "simpletrim".to_string(),
        "long".into(),
        "-i".into(),
        args.r1.display().to_string(),
        "-o".into(),
        trimmed_name(&args.r1),
    ];
    for s in &suggestions {
        cmd.extend(s.flags.iter().cloned());
    }
    let command = cmd
        .iter()
        .map(|a| shell_quote(a))
        .collect::<Vec<_>>()
        .join(" ");

    let profile = |bins: &[(f64, u64); PROFILE_BINS]| -> Vec<Value> {
        bins.iter().map(|&(e, n)| q_or_null(e, n)).collect()
    };
    let yields: Vec<Value> = YIELD_LENGTHS
        .iter()
        .map(|&l| {
            let b = rs.lengths.bases_at_least(l);
            json!({ "min_length": l, "bases": b, "fraction": round(b as f64 / rs.base.bases.max(1) as f64, 5) })
        })
        .collect();
    let mut reads_json = rs.to_json();
    reads_json["yield_by_min_length"] = json!(yields);
    reads_json["quality_profile"] = json!({
        "averaging": "error_probability",
        "bin_size": PROFILE_BIN,
        "from_start": profile(&st.from_start),
        "from_end": profile(&st.from_end),
    });
    reads_json["poly_x_tails"] = json!({
        "min_len": POLY_X_MIN_LEN,
        "reads": poly_x,
        "fraction": round(rate(poly_x), 5),
        "by_base": { "A": st.poly_x[0], "C": st.poly_x[1], "G": st.poly_x[2], "T": st.poly_x[3] },
    });
    reads_json["internal_low_quality"] = json!({
        "window": LOW_QUALITY_WINDOW,
        "threshold": LOW_QUALITY_THRESHOLD,
        "reads": st.internal_low_quality,
        "fraction": round(rate(st.internal_low_quality), 5),
        "bases": st.internal_low_quality_bases,
    });
    reads_json["rna"] = json!(rna);
    let end_json = |e: &crate::long_detect::EndDetection, fraction: f64, found: bool| {
        json!({
            "detected": e.adapter,
            "too_short_candidate": e.too_short,
            "reads_with_adapter_fraction": found.then(|| round(fraction, 5)),
        })
    };
    let mut platform = header_json(header);
    platform["run"] = run_info(first_name);
    let report = json!({
        "version": env!("CARGO_PKG_VERSION"),
        "status": "ok",
        "mode": "long",
        "inputs": { "r1": file_json(&args.r1), "r2": null },
        "platform": platform,
        "sample": {
            "records": total,
            "bases": rs.base.bases,
            "whole_file": c.whole_file,
            "estimated_total_records": c.estimated_total,
        },
        "reads": reads_json,
        "adapters": {
            "detection_reads": d.reads,
            "check_sample_reads": c.sample.len(),
            "start": end_json(&d.start, at_start, start.is_some()),
            "end": end_json(&d.end, at_end, end.is_some()),
            "inside_reads_fraction": any_adapter.then(|| round(inside, 5)),
        },
        "warnings": warnings,
        "notes": notes,
        "suggestions": suggestions.iter().map(|s| json!({"flags": s.flags, "reason": s.reason})).collect::<Vec<_>>(),
        "suggested_command": command,
    });
    let text = render_text(args, &report);
    emit(log, &text, &report)
}

fn pct(v: &Value) -> String {
    v.as_f64()
        .map_or("-".into(), |f| format!("{:.2}%", f * 100.0))
}

fn num(v: &Value) -> String {
    match v {
        Value::Number(n) if n.is_u64() => group(n.as_u64().unwrap_or(0)),
        Value::Number(n) => format!("{:.1}", n.as_f64().unwrap_or(0.0)),
        Value::Null => "-".into(),
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

fn qual(v: &Value) -> String {
    v.as_f64().map_or("-".into(), |q| format!("Q{q:.1}"))
}

fn render_text(args: &AnalyzeArgs, r: &Value) -> String {
    let mut t = String::new();
    let mut line = |s: String| {
        t.push_str(&s);
        t.push('\n');
    };
    let sample = &r["sample"];
    let scope = if sample["whole_file"].as_bool() == Some(true) {
        "whole file".to_string()
    } else {
        match sample["estimated_total_records"].as_u64() {
            Some(est) => format!("first part of ~{} reads", group(est)),
            None => "first part of the file".into(),
        }
    };
    line(format!(
        "simpletrim analyze: {} long reads inspected ({scope})",
        num(&sample["records"])
    ));
    line(String::new());
    line("Input".into());
    let f = &r["inputs"]["r1"];
    let size = f["file_size"]
        .as_u64()
        .map(|b| format!(", {:.2} GB", b as f64 / 1e9))
        .unwrap_or_default();
    line(format!(
        "  File: {} ({}{size})",
        args.r1.display(),
        if f["compressed"].as_bool() == Some(true) {
            "gzip"
        } else {
            "plain"
        }
    ));
    let p = &r["platform"];
    let mut platform = p["name"].as_str().unwrap_or("unknown").to_string();
    if let Some(run) = p["run"].as_object() {
        for (key, label) in [
            ("basecall_model", "basecall model"),
            ("flow_cell", "flow cell"),
            ("barcode", "barcode"),
            ("run_id", "run"),
        ] {
            if let Some(v) = run.get(key).and_then(Value::as_str) {
                platform.push_str(&format!("; {label} {v}"));
            }
        }
    }
    line(format!("  Platform: {platform}"));
    line(String::new());

    let rd = &r["reads"];
    let row = |label: &str, value: String| format!("  {label:<38}{value}");
    line("Reads".into());
    line(row(
        "reads / bases",
        format!("{} / {}", num(&rd["reads"]), num(&rd["bases"])),
    ));
    let len = &rd["length"];
    line(row(
        "length min / median / max",
        format!(
            "{} / {} / {}",
            num(&len["min"]),
            num(&len["median"]),
            num(&len["max"])
        ),
    ));
    line(row(
        "mean length / N50 / N90",
        format!(
            "{} / {} / {}",
            num(&len["mean"]),
            num(&len["n50"]),
            num(&len["n90"])
        ),
    ));
    let yields: Vec<String> = rd["yield_by_min_length"]
        .as_array()
        .map(|a| a.iter().map(|y| pct(&y["fraction"])).collect())
        .unwrap_or_default();
    line(row(
        "bases in reads >= 1 / 5 / 10 / 20 kb",
        yields.join(" / "),
    ));
    line(row(
        "mean quality (error probability)",
        qual(&rd["mean_quality"]),
    ));
    line(row(
        "median read quality",
        rd["median_read_quality"]
            .as_u64()
            .map_or("-".into(), |q| format!("Q{q}")),
    ));
    let total = rd["reads"].as_f64().unwrap_or(0.0).max(1.0);
    let levels: Vec<String> = ["q7", "q10", "q15", "q20"]
        .iter()
        .map(|k| {
            let n = rd["reads_with_mean_quality_at_least"][k]
                .as_f64()
                .unwrap_or(0.0);
            format!("{:.1}%", n * 100.0 / total)
        })
        .collect();
    line(row("reads >= Q7 / Q10 / Q15 / Q20", levels.join(" / ")));
    line(row(
        "Q20 / Q30 bases",
        format!("{} / {}", pct(&rd["q20_rate"]), pct(&rd["q30_rate"])),
    ));
    line(row("GC content", pct(&rd["gc_content"])));
    line(row(
        "poly-X tails (>= 10 bp)",
        pct(&rd["poly_x_tails"]["fraction"]),
    ));
    line(row(
        "low-quality stretches inside reads",
        pct(&rd["internal_low_quality"]["fraction"]),
    ));
    line(String::new());

    let prof = &rd["quality_profile"];
    line(format!(
        "Quality along the reads (error-probability mean per {} bases)",
        num(&prof["bin_size"])
    ));
    for (label, key) in [("from start", "from_start"), ("from end", "from_end")] {
        let bins: Vec<String> = prof[key]
            .as_array()
            .map(|a| a.iter().map(qual).collect())
            .unwrap_or_default();
        line(format!("  {label:<12}{}", bins.join(" ")));
    }
    line(String::new());

    let ad = &r["adapters"];
    line(format!(
        "Adapters (detected on the first {} reads)",
        num(&ad["detection_reads"])
    ));
    for which in ["start", "end"] {
        let a = &ad[which];
        let desc = match a["detected"].as_str() {
            Some(s) => format!(
                "{s} - in {} of {} checked reads",
                pct(&a["reads_with_adapter_fraction"]),
                num(&ad["check_sample_reads"])
            ),
            None => "none detected".into(),
        };
        line(format!("  {which:<7}{desc}"));
    }
    if !ad["inside_reads_fraction"].is_null() {
        line(format!(
            "  inside {} of checked reads",
            pct(&ad["inside_reads_fraction"])
        ));
    }
    line(String::new());

    for (title, key) in [("Warnings", "warnings"), ("Notes", "notes")] {
        if let Some(items) = r[key].as_array().filter(|a| !a.is_empty()) {
            line(title.into());
            for i in items {
                line(format!("  - {}", i.as_str().unwrap_or("")));
            }
            line(String::new());
        }
    }
    line("Suggestions".into());
    match r["suggestions"].as_array().filter(|a| !a.is_empty()) {
        Some(items) => {
            for s in items {
                let flags: Vec<&str> = s["flags"]
                    .as_array()
                    .map(|f| f.iter().filter_map(|x| x.as_str()).collect())
                    .unwrap_or_default();
                line(format!(
                    "  {}: {}",
                    flags.join(" "),
                    s["reason"].as_str().unwrap_or("")
                ));
            }
        }
        None => line("  the defaults look fine".into()),
    }
    line(String::new());
    line("Suggested command".into());
    line(format!(
        "  {}",
        r["suggested_command"].as_str().unwrap_or("")
    ));
    t
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn dorado_and_minknow_headers() {
        let dorado = b"@3e50906d-90df-4f4c-8087-f70d7a3b4b5a\tqs:f:14.1\tch:i:91\tRG:Z:96dcd2f3-53f7-4a88-a69c-ebaa38c687cb_dna_r10.4.1_e8.2_400bps_hac@v5.2.0_barcode01\tPU:Z:FBH77215\tal:Z:barcode01";
        let v = run_info(dorado);
        assert_eq!(v["run_id"], "96dcd2f3-53f7-4a88-a69c-ebaa38c687cb");
        assert_eq!(v["basecall_model"], "dna_r10.4.1_e8.2_400bps_hac@v5.2.0");
        assert_eq!(v["flow_cell"], "FBH77215");
        assert_eq!(v["barcode"], "barcode01");
        let minknow = b"@0a1b2c3d-1234-4abc-9def-0123456789ab runid=abc123 read=1 ch=5 flow_cell_id=FAQ12345 basecall_model_version_id=dna_r9.4.1_450bps_hac barcode=barcode07";
        let v = run_info(minknow);
        assert_eq!(v["run_id"], "abc123");
        assert_eq!(v["flow_cell"], "FAQ12345");
        assert_eq!(v["basecall_model"], "dna_r9.4.1_450bps_hac");
        assert_eq!(v["barcode"], "barcode07");
    }
}
