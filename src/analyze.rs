//! The `analyze` subcommand: inspects the inputs and suggests a `simpletrim short` command,
//! or for long reads a `simpletrim long` command (see [`long_reads`]).
//!
//! The report goes to stdout (text, or JSON with `--json`); progress messages go to stderr.

mod long_reads;

use std::io::Write;
use std::path::Path;
use std::sync::atomic::Ordering;
use std::sync::mpsc::{Receiver, sync_channel};

use rayon::prelude::*;
use serde_json::{Value, json};

use crate::Result;
use crate::adapter::{MatchParams, trim_single};
use crate::cli::AnalyzeArgs;
use crate::detect::{self, Detection};
use crate::fastq::{Chunk, Record, feed_chunks, mate_names_match};
use crate::io::{self, is_gz};
use crate::log::Logger;
use crate::overlap::{self, OverlapParams};
use crate::platform::{HeaderInfo, Platform, parse_header};
use crate::polyg::poly_g_trim_len;
use crate::stats::{group, round};

const CHUNK_RECORDS: usize = 8192;
/// Insert sizes are counted up to this value.
const MAX_INSERT: usize = 1000;
/// Poly-G tails at least this long are counted (fastp default).
const POLY_G_MIN_LEN: usize = 10;
/// Motifs counted by exact occurrence (in the style of FastQC's adapter content).
const MOTIFS: [(&str, &str, &str); 5] = [
    (
        "illumina_universal",
        "Illumina universal adapter",
        "AGATCGGAAGAG",
    ),
    ("nextera", "Nextera / Tn5 adapter", "CTGTCTCTTATA"),
    (
        "small_rna_3p",
        "Illumina small RNA 3' adapter",
        "TGGAATTCTCGG",
    ),
    ("mgi", "MGI / DNBSEQ adapter", "AAGTCGGAGGCC"),
    ("poly_a", "poly-A", "AAAAAAAAAAAA"),
];
/// A suggestion is made if a feature is seen in at least this fraction of reads.
const POLY_G_SUGGEST_RATE: f64 = 0.001;
const ADAPTER_SUGGEST_RATE: f64 = 0.01;

#[derive(Clone)]
struct MateStats {
    reads: u64,
    bases: u64,
    len_hist: Vec<u64>,
    qual_hist: Vec<u64>,
    qual_sum: u64,
    q20: u64,
    q30: u64,
    pos_qual: Vec<u64>,
    pos_count: Vec<u64>,
    /// Per position counts of A, C, G, T, other.
    pos_comp: Vec<[u64; 5]>,
    poly_g_reads: u64,
    motif_reads: [u64; MOTIFS.len()],
}

impl MateStats {
    fn new() -> Self {
        MateStats {
            reads: 0,
            bases: 0,
            len_hist: Vec::new(),
            qual_hist: vec![0; 256],
            qual_sum: 0,
            q20: 0,
            q30: 0,
            pos_qual: Vec::new(),
            pos_count: Vec::new(),
            pos_comp: Vec::new(),
            poly_g_reads: 0,
            motif_reads: [0; MOTIFS.len()],
        }
    }

    fn add(&mut self, rec: &Record) {
        let len = rec.seq.len();
        self.reads += 1;
        self.bases += len as u64;
        if self.len_hist.len() <= len {
            self.len_hist.resize(len + 1, 0);
        }
        self.len_hist[len] += 1;
        if self.pos_count.len() < len {
            self.pos_qual.resize(len, 0);
            self.pos_count.resize(len, 0);
            self.pos_comp.resize(len, [0; 5]);
        }
        for (i, (&b, &q)) in rec.seq.iter().zip(rec.qual).enumerate() {
            self.qual_hist[q as usize] += 1;
            let phred = q.saturating_sub(33) as u64;
            self.qual_sum += phred;
            self.q20 += (phred >= 20) as u64;
            self.q30 += (phred >= 30) as u64;
            self.pos_qual[i] += phred;
            self.pos_count[i] += 1;
            let base = match b {
                b'A' | b'a' => 0,
                b'C' | b'c' => 1,
                b'G' | b'g' => 2,
                b'T' | b't' => 3,
                _ => 4,
            };
            self.pos_comp[i][base] += 1;
        }
        if poly_g_trim_len(rec.seq, POLY_G_MIN_LEN).is_some() {
            self.poly_g_reads += 1;
        }
        for (k, (_, _, motif)) in MOTIFS.iter().enumerate() {
            let m = motif.as_bytes();
            if rec.seq.len() >= m.len() && rec.seq.windows(m.len()).any(|w| w == m) {
                self.motif_reads[k] += 1;
            }
        }
    }

    fn merge(&mut self, o: &MateStats) {
        fn add_vec<T: Copy + Default>(a: &mut Vec<T>, b: &[T], f: impl Fn(&mut T, &T)) {
            if a.len() < b.len() {
                a.resize(b.len(), T::default());
            }
            for (x, y) in a.iter_mut().zip(b) {
                f(x, y);
            }
        }
        self.reads += o.reads;
        self.bases += o.bases;
        self.qual_sum += o.qual_sum;
        self.q20 += o.q20;
        self.q30 += o.q30;
        self.poly_g_reads += o.poly_g_reads;
        add_vec(&mut self.len_hist, &o.len_hist, |x, y| *x += y);
        add_vec(&mut self.qual_hist, &o.qual_hist, |x, y| *x += y);
        add_vec(&mut self.pos_qual, &o.pos_qual, |x, y| *x += y);
        add_vec(&mut self.pos_count, &o.pos_count, |x, y| *x += y);
        add_vec(&mut self.pos_comp, &o.pos_comp, |x, y| {
            for i in 0..5 {
                x[i] += y[i];
            }
        });
        for (x, y) in self.motif_reads.iter_mut().zip(&o.motif_reads) {
            *x += y;
        }
    }

    fn rate(&self, n: u64) -> f64 {
        if self.reads == 0 {
            0.0
        } else {
            n as f64 / self.reads as f64
        }
    }

    fn base_rate(&self, n: u64) -> f64 {
        if self.bases == 0 {
            0.0
        } else {
            n as f64 / self.bases as f64
        }
    }

    fn composition(&self) -> [u64; 5] {
        let mut total = [0u64; 5];
        for p in &self.pos_comp {
            for i in 0..5 {
                total[i] += p[i];
            }
        }
        total
    }

    /// Mean phred quality per position, for positions covered by at least half of the reads.
    fn cycle_means(&self) -> Vec<f64> {
        self.pos_qual
            .iter()
            .zip(&self.pos_count)
            .take_while(|&(_, &n)| n * 2 >= self.reads && n > 0)
            .map(|(&q, &n)| q as f64 / n as f64)
            .collect()
    }

    /// 1-based positions among the first 20 where A/T or G/C frequencies differ by > 10 points.
    fn biased_positions(&self) -> Vec<usize> {
        self.pos_comp
            .iter()
            .take(20)
            .enumerate()
            .filter_map(|(i, c)| {
                let n = (c[0] + c[1] + c[2] + c[3]) as f64;
                if n == 0.0 {
                    return None;
                }
                let f = |k: usize| c[k] as f64 / n;
                ((f(0) - f(3)).abs() > 0.1 || (f(1) - f(2)).abs() > 0.1).then_some(i + 1)
            })
            .collect()
    }
}

#[derive(Clone, Default)]
struct PairStats {
    pairs: u64,
    name_mismatches: u64,
    first_mismatch: Option<(u64, String, String)>,
    overlapped: u64,
    read_through: u64,
    insert_hist: Vec<u64>,
}

impl PairStats {
    fn merge(&mut self, o: &PairStats) {
        self.pairs += o.pairs;
        self.name_mismatches += o.name_mismatches;
        if self.first_mismatch.is_none() {
            self.first_mismatch.clone_from(&o.first_mismatch);
        }
        self.overlapped += o.overlapped;
        self.read_through += o.read_through;
        if self.insert_hist.len() < o.insert_hist.len() {
            self.insert_hist.resize(o.insert_hist.len(), 0);
        }
        for (x, y) in self.insert_hist.iter_mut().zip(&o.insert_hist) {
            *x += y;
        }
    }
}

#[derive(Clone)]
struct Partial {
    mates: [MateStats; 2],
    pair: PairStats,
}

impl Partial {
    fn new() -> Self {
        Partial {
            mates: [MateStats::new(), MateStats::new()],
            pair: PairStats::default(),
        }
    }

    fn merge(&mut self, o: &Partial) {
        self.mates[0].merge(&o.mates[0]);
        self.mates[1].merge(&o.mates[1]);
        self.pair.merge(&o.pair);
    }
}

fn analyze_chunk(c1: &Chunk, c2: Option<&Chunk>) -> Partial {
    let mut p = Partial::new();
    let params = OverlapParams::default();
    let mut rc_buf = Vec::new();
    for i in 0..c1.len() {
        let r1 = c1.get(i);
        p.mates[0].add(&r1);
        let Some(c2) = c2 else { continue };
        let r2 = c2.get(i);
        p.mates[1].add(&r2);
        let ps = &mut p.pair;
        ps.pairs += 1;
        if !mate_names_match(&r1, &r2) {
            ps.name_mismatches += 1;
            if ps.first_mismatch.is_none() {
                ps.first_mismatch = Some((
                    c1.first_record + i as u64 + 1,
                    String::from_utf8_lossy(r1.id()).into_owned(),
                    String::from_utf8_lossy(r2.id()).into_owned(),
                ));
            }
        }
        let ov = overlap::analyze(r1.seq, r2.seq, &params, &mut rc_buf);
        if let Some(size) = overlap::insert_size(r1.seq.len(), r2.seq.len(), &ov) {
            ps.overlapped += 1;
            ps.read_through += (ov.offset < 0) as u64;
            let size = size.min(MAX_INSERT);
            if ps.insert_hist.len() <= size {
                ps.insert_hist.resize(size + 1, 0);
            }
            ps.insert_hist[size] += 1;
        }
    }
    p
}

/// Value at quantile `q` of a histogram (index = value).
fn hist_quantile(hist: &[u64], q: f64) -> Option<usize> {
    let total: u64 = hist.iter().sum();
    if total == 0 {
        return None;
    }
    let target = (q * total as f64).ceil().max(1.0) as u64;
    let mut acc = 0;
    for (v, &n) in hist.iter().enumerate() {
        acc += n;
        if acc >= target {
            return Some(v);
        }
    }
    None
}

fn hist_mode(hist: &[u64]) -> Option<usize> {
    let (v, &n) = hist
        .iter()
        .enumerate()
        .max_by_key(|&(v, &n)| (n, std::cmp::Reverse(v)))?;
    (n > 0).then_some(v)
}

fn hist_min_max(hist: &[u64]) -> Option<(usize, usize)> {
    let min = hist.iter().position(|&n| n > 0)?;
    let max = hist.iter().rposition(|&n| n > 0)?;
    Some((min, max))
}

fn quality_encoding(qual_hist: &[u64]) -> (&'static str, Option<u8>, Option<u8>) {
    let min = qual_hist.iter().position(|&n| n > 0).map(|c| c as u8);
    let max = qual_hist.iter().rposition(|&n| n > 0).map(|c| c as u8);
    let encoding = match min {
        None => "unknown",
        Some(c) if c < 59 => "phred+33",
        Some(c) if c >= 64 => "phred+64 (likely)",
        Some(_) => "phred+33 (likely)",
    };
    (encoding, min, max)
}

/// Output name for a suggested command: `x.fastq.gz` becomes `x.trimmed.fastq.gz`.
/// The directory part is kept exactly as given, separators included.
fn trimmed_name(p: &Path) -> String {
    let full = p.to_string_lossy();
    let name = p
        .file_name()
        .map(|n| n.to_string_lossy().into_owned())
        .unwrap_or_default();
    let dir = full.strip_suffix(name.as_str()).unwrap_or("");
    let lower = name.to_ascii_lowercase();
    for ext in [".fastq.gz", ".fq.gz", ".fastq", ".fq"] {
        if lower.ends_with(ext) && lower.len() > ext.len() {
            let (stem, ext) = name.split_at(name.len() - ext.len());
            return format!("{dir}{stem}.trimmed{ext}");
        }
    }
    format!("{dir}{name}.trimmed.fastq.gz")
}

fn shell_quote(s: &str) -> String {
    if !s.is_empty()
        && s.chars()
            .all(|c| c.is_ascii_alphanumeric() || "-_./\\:=+,@%".contains(c))
    {
        s.to_string()
    } else {
        format!("\"{}\"", s.replace('"', "\\\""))
    }
}

fn file_json(path: &Path) -> Value {
    json!({
        "path": path.display().to_string(),
        "compressed": is_gz(path),
        "file_size": std::fs::metadata(path).map(|m| m.len()).ok(),
    })
}

fn header_json(h: &HeaderInfo) -> Value {
    json!({
        "platform": h.platform,
        "name": h.platform.label(),
        "header_format": h.format,
        "instrument_id": h.instrument_id,
        "instrument_model": h.instrument_model,
        "two_color": h.two_color,
        "read_number": h.read_number,
    })
}

fn peek_header(path: &Path, label: &str) -> Result<Option<HeaderInfo>> {
    let mut input = io::open_fastq(path, label)?;
    Ok(input
        .reader
        .read_chunk(1, usize::MAX)?
        .map(|c| parse_header(c.get(0).name)))
}

fn emit(log: &Logger, text: &str, json: &Value) -> Result<()> {
    let mut out = std::io::stdout().lock();
    let result = if log.is_json() {
        writeln!(
            out,
            "{}",
            serde_json::to_string_pretty(json).expect("serializable")
        )
    } else {
        write!(out, "{text}")
    };
    result
        .and_then(|_| out.flush())
        .map_err(|e| format!("cannot write report: {e}"))
}

/// Reads that are this long (median of the first reads) are analysed as long reads even if
/// the header does not tell the platform.
const LONG_READ_MEDIAN: usize = 1000;
const PEEK_READS: usize = 100;

/// Header of the first read, its name, and the median length of the first reads.
fn peek_start(path: &Path, label: &str) -> Result<Option<(HeaderInfo, Vec<u8>, usize)>> {
    let mut input = io::open_fastq(path, label)?;
    let Some(chunk) = input.reader.read_chunk(PEEK_READS, usize::MAX)? else {
        return Ok(None);
    };
    let mut lens: Vec<usize> = chunk.iter().map(|r| r.seq.len()).collect();
    lens.sort_unstable();
    let name = chunk.get(0).name.to_vec();
    Ok(Some((parse_header(&name), name, lens[lens.len() / 2])))
}

fn receive(rx: Option<&Receiver<Result<Chunk>>>) -> Result<Option<Chunk>> {
    match rx.map(|rx| rx.recv()) {
        Some(Ok(chunk)) => chunk.map(Some),
        _ => Ok(None),
    }
}

struct Collected {
    stats: Partial,
    samples: [Vec<Vec<u8>>; 2],
    /// Set if R1 and R2 ended at different read counts.
    count_mismatch: Option<String>,
    whole_file: bool,
    estimated_total: Option<u64>,
}

fn collect(args: &AnalyzeArgs, threads: usize, log: &Logger) -> Result<Collected> {
    let in1 = io::open_fastq(&args.r1, "R1")?;
    let in2 = args
        .r2
        .as_deref()
        .map(|p| io::open_fastq(p, "R2"))
        .transpose()?;
    let bytes_read = in1.bytes_read.clone();
    let file_size = in1.file_size;
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| format!("cannot start worker threads: {e}"))?;
    let batch_size = threads * 2;

    std::thread::scope(|scope| {
        let (tx1, rx1) = sync_channel(batch_size);
        scope
            .spawn(move || feed_chunks(in1.reader, tx1, args.max_reads, CHUNK_RECORDS, usize::MAX));
        let rx2 = in2.map(|in2| {
            let (tx2, rx2) = sync_channel(batch_size);
            scope.spawn(move || {
                feed_chunks(in2.reader, tx2, args.max_reads, CHUNK_RECORDS, usize::MAX)
            });
            rx2
        });
        let mut stats = Partial::new();
        let mut samples: [Vec<Vec<u8>>; 2] = [Vec::new(), Vec::new()];
        let mut sample_bases = [0usize; 2];
        let mut count_mismatch = None;
        let mut next_report = 1_000_000u64;
        'outer: loop {
            let mut batch = Vec::with_capacity(batch_size);
            let mut ended = false;
            while batch.len() < batch_size {
                let c1 = receive(Some(&rx1))?;
                let c2 = receive(rx2.as_ref())?;
                match (c1, c2) {
                    (None, None) => {
                        ended = true;
                        break;
                    }
                    (Some(c1), None) if rx2.is_none() => batch.push((c1, None)),
                    (Some(c1), Some(c2)) if c1.len() == c2.len() => batch.push((c1, Some(c2))),
                    (c1, c2) => {
                        let end =
                            |c: &Option<Chunk>| c.as_ref().map(|c| c.first_record + c.len() as u64);
                        count_mismatch = Some(format!(
                            "R1 and R2 have different numbers of reads (R1 has {} or more, R2 {} or more at the point where they diverge)",
                            group(end(&c1).unwrap_or(0)),
                            group(end(&c2).unwrap_or(0)),
                        ));
                        ended = true;
                        break;
                    }
                }
            }
            for (c1, c2) in &batch {
                for (m, chunk) in [Some(c1), c2.as_ref()].into_iter().enumerate() {
                    let Some(chunk) = chunk else { continue };
                    for rec in chunk.iter() {
                        if samples[m].len() >= detect::SAMPLE_READS
                            || sample_bases[m] >= detect::SAMPLE_BASES
                        {
                            break;
                        }
                        sample_bases[m] += rec.seq.len();
                        samples[m].push(rec.seq.to_vec());
                    }
                }
            }
            let partials: Vec<Partial> = pool.install(|| {
                batch
                    .par_iter()
                    .map(|(c1, c2)| analyze_chunk(c1, c2.as_ref()))
                    .collect()
            });
            for p in &partials {
                stats.merge(p);
            }
            if stats.mates[0].reads >= next_report {
                log.progress(
                    &format!(
                        "read {} {}",
                        group(stats.mates[0].reads),
                        if rx2.is_some() { "pairs" } else { "reads" }
                    ),
                    json!({ "records": stats.mates[0].reads }),
                );
                next_report += 1_000_000;
            }
            if ended {
                break 'outer;
            }
        }
        drop(rx1);
        drop(rx2);
        let reads = stats.mates[0].reads;
        let consumed = bytes_read.load(Ordering::Relaxed);
        let whole_file = args.max_reads == 0
            || reads < args.max_reads
            || (file_size > 0 && consumed >= file_size);
        let estimated_total = (!whole_file && consumed > 0)
            .then(|| (reads as f64 * file_size as f64 / consumed as f64).round() as u64);
        Ok(Collected {
            stats,
            samples,
            count_mismatch,
            whole_file,
            estimated_total,
        })
    })
}

struct Suggestion {
    flags: Vec<String>,
    reason: String,
}

/// Runs the `analyze` subcommand.
pub fn run(args: &AnalyzeArgs, threads: usize, log: &Logger) -> Result<()> {
    io::check_input(&args.r1, "R1 input")?;
    let paired = args.r2.is_some();
    if let Some(r2) = &args.r2 {
        io::check_input(r2, "R2 input")?;
        if io::same_path(&args.r1, r2) {
            return Err("--r1 and --r2 point to the same file".into());
        }
    }
    let Some((h1, first_name, median_len)) = peek_start(&args.r1, "R1")? else {
        return Err(format!("{} contains no reads", args.r1.display()));
    };
    if h1.platform.is_long_read()
        || (h1.platform == Platform::Unknown && median_len >= LONG_READ_MEDIAN)
    {
        return long_reads::run(args, &h1, &first_name, threads, log);
    }
    let h2 = match &args.r2 {
        Some(p) => peek_header(p, "R2")?,
        None => None,
    };

    let unit = if paired { "pairs" } else { "reads" };
    log.info(&format!(
        "Reading {} ...",
        if args.max_reads == 0 {
            format!("all {unit}")
        } else {
            format!("up to {} {unit}", group(args.max_reads))
        }
    ));
    let c = collect(args, threads, log)?;
    let mates = if paired { 2 } else { 1 };

    log.info("Detecting adapters ...");
    let detections: Vec<Detection> = std::thread::scope(|s| {
        let handles: Vec<_> = c.samples[..mates]
            .iter()
            .map(|sample| s.spawn(move || detect::detect(sample, 0)))
            .collect();
        handles
            .into_iter()
            .map(|h| h.join().expect("adapter detection panicked"))
            .collect()
    });
    // fraction of sampled reads in which the detected adapter is found
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(threads)
        .build()
        .map_err(|e| format!("cannot start worker threads: {e}"))?;
    let prevalence: Vec<Option<f64>> = (0..mates)
        .map(|m| {
            let adapter = detections[m].adapter.as_ref()?.sequence.as_bytes().to_vec();
            let sample = &c.samples[m];
            if sample.is_empty() {
                return None;
            }
            let params = MatchParams::default();
            let hits = pool.install(|| {
                sample
                    .par_iter()
                    .filter(|s| trim_single(s, &adapter, &params).is_some_and(|n| n + 4 < s.len()))
                    .count()
            });
            Some(hits as f64 / sample.len() as f64)
        })
        .collect();

    let st = &c.stats;
    let total = st.mates[0].reads;
    let mut warnings: Vec<String> = Vec::new();
    let mut notes: Vec<String> = Vec::new();
    let mut suggestions: Vec<Suggestion> = Vec::new();

    if total == 0 {
        return Err(format!("{} contains no reads", args.r1.display()));
    }

    // pairing
    if let Some(msg) = &c.count_mismatch {
        warnings.push(format!("{msg}; 'simpletrim short' will stop with an error"));
    }
    let ps = &st.pair;
    if paired && ps.name_mismatches > 0 {
        let (n, a, b) = ps.first_mismatch.clone().unwrap_or_default();
        warnings.push(format!(
            "read names differ in {} of {} pairs (first at pair {n}: '{a}' vs '{b}'); \
             'simpletrim short' aborts on this unless --no-name-check is given",
            group(ps.name_mismatches),
            group(ps.pairs)
        ));
    }
    if let (Some(n1), Some(n2)) = (h1.read_number, h2.as_ref().and_then(|h| h.read_number)) {
        if n1 == 2 && n2 == 1 {
            warnings.push(
                "headers mark --r1 as read 2 and --r2 as read 1; the files may be swapped".into(),
            );
        } else if n1 == n2 {
            warnings.push(format!(
                "both files are marked as read {n1} in their headers"
            ));
        }
    }
    if let Some(h2) = &h2
        && (h2.platform != h1.platform || h2.instrument_id != h1.instrument_id)
    {
        warnings.push("R1 and R2 headers suggest different instruments".into());
    }

    // per-mate figures
    let mut mate_json = Vec::new();
    let mut poly_g_rate_max: f64 = 0.0;
    let mut low_tail_cycle: Option<(usize, usize)> = None;
    let mut low_front = false;
    for (m, ms) in st.mates[..mates].iter().enumerate() {
        let (min_len, max_len) = hist_min_max(&ms.len_hist).unwrap_or((0, 0));
        let median_len = hist_quantile(&ms.len_hist, 0.5).unwrap_or(0);
        let (encoding, qmin, qmax) = quality_encoding(&ms.qual_hist);
        if encoding.starts_with("phred+64") {
            warnings.push(format!(
                "R{} qualities look like phred+64, which simpletrim does not support",
                m + 1
            ));
        }
        let distinct = ms.qual_hist.iter().filter(|&&n| n > 0).count();
        let comp = ms.composition();
        let cycles = ms.cycle_means();
        let below = |q: f64| cycles.iter().position(|&x| x < q).map(|i| i + 1);
        let tail_start = cycles.len().saturating_sub(cycles.len() / 3).max(5);
        if let Some(i) = cycles.iter().skip(tail_start).position(|&x| x < 20.0) {
            let cycle = tail_start + i + 1;
            low_tail_cycle = Some(
                low_tail_cycle.map_or((m, cycle), |p| if cycle < p.1 { (m, cycle) } else { p }),
            );
        }
        if cycles.iter().take(3).any(|&x| x < 20.0) {
            low_front = true;
        }
        let poly_g_rate = ms.rate(ms.poly_g_reads);
        poly_g_rate_max = poly_g_rate_max.max(poly_g_rate);
        let motifs: serde_json::Map<String, Value> = MOTIFS
            .iter()
            .enumerate()
            .map(|(k, (key, name, seq))| {
                (
                    key.to_string(),
                    json!({ "name": name, "motif": seq, "reads": ms.motif_reads[k], "fraction": round(ms.rate(ms.motif_reads[k]), 5) }),
                )
            })
            .collect();
        let d = &detections[m];
        mate_json.push(json!({
            "reads": ms.reads,
            "bases": ms.bases,
            "length": {
                "min": min_len,
                "max": max_len,
                "mean": round(ms.bases as f64 / ms.reads.max(1) as f64, 2),
                "median": median_len,
                "mode": hist_mode(&ms.len_hist),
                "fixed": min_len == max_len,
            },
            "quality": {
                "encoding": encoding,
                "min_char": qmin.map(|c| (c as char).to_string()),
                "max_char": qmax.map(|c| (c as char).to_string()),
                "distinct_values": distinct,
                "binned": distinct <= 8,
                "mean": round(ms.base_rate(ms.qual_sum), 2),
                "q20_rate": round(ms.base_rate(ms.q20), 5),
                "q30_rate": round(ms.base_rate(ms.q30), 5),
                "first_cycle_below_q30": below(30.0),
                "first_cycle_below_q20": below(20.0),
                "per_cycle_mean": cycles.iter().map(|&q| round(q, 2)).collect::<Vec<_>>(),
            },
            "gc_content": round(ms.base_rate(comp[1] + comp[2]), 5),
            "n_rate": round(ms.base_rate(comp[4]), 6),
            "poly_g": { "min_len": POLY_G_MIN_LEN, "reads": ms.poly_g_reads, "fraction": round(poly_g_rate, 5) },
            "motifs": motifs,
            "composition_bias_cycles": ms.biased_positions(),
            "adapter": {
                "detected": d.adapter.as_ref().map(|a| json!({ "sequence": a.sequence, "name": a.name })),
                "unconfirmed_candidate": d.candidate,
                "sample_reads": c.samples[m].len(),
                "reads_with_adapter_fraction": prevalence[m].map(|f| round(f, 5)),
            },
        }));
    }

    // poly-G
    if poly_g_rate_max >= POLY_G_SUGGEST_RATE {
        if h1.two_color {
            notes.push(format!(
                "poly-G tails in up to {:.2}% of reads; poly-G trimming is enabled automatically for this two-colour instrument",
                poly_g_rate_max * 100.0
            ));
        } else {
            suggestions.push(Suggestion {
                flags: vec!["--trim-poly-g".into()],
                reason: format!(
                    "poly-G tails in up to {:.2}% of reads (not enabled automatically for this instrument)",
                    poly_g_rate_max * 100.0
                ),
            });
        }
    } else if h1.two_color {
        notes.push("two-colour instrument: poly-G trimming will be enabled automatically, although few poly-G tails were seen".into());
    }

    // adapters
    let adapter_seq = |m: usize| detections[m].adapter.as_ref().map(|a| a.sequence.clone());
    let frequent = |m: usize| prevalence[m].is_some_and(|f| f >= ADAPTER_SUGGEST_RATE);
    if paired {
        let read_through = if ps.pairs > 0 {
            ps.read_through as f64 / ps.pairs as f64
        } else {
            0.0
        };
        if read_through > 0.0 {
            notes.push(format!(
                "{:.2}% of pairs overlap with adapter read-through; these are trimmed by overlap analysis automatically",
                read_through * 100.0
            ));
        }
        if (frequent(0) || frequent(1)) && (adapter_seq(0).is_some() || adapter_seq(1).is_some()) {
            let mut flags = Vec::new();
            for m in 0..2 {
                if let Some(a) = adapter_seq(m) {
                    flags.push(format!("--adapter-r{}", m + 1));
                    flags.push(a);
                }
            }
            suggestions.push(Suggestion {
                flags,
                reason:
                    "detected adapters, used for pairs whose reads cannot be aligned to each other"
                        .into(),
            });
        }
    } else if let Some(a) = adapter_seq(0) {
        suggestions.push(Suggestion {
            flags: vec!["--adapter-r1".into(), a],
            reason: "detected adapter (same as the default auto-detection, but skips the detection pass)".into(),
        });
    } else {
        let motif = |k: usize| st.mates[0].rate(st.mates[0].motif_reads[k]);
        if motif(0) >= ADAPTER_SUGGEST_RATE {
            suggestions.push(Suggestion {
                flags: vec!["--adapter-preset".into(), "truseq".into()],
                reason: format!(
                    "Illumina universal adapter motif in {:.2}% of reads",
                    motif(0) * 100.0
                ),
            });
        } else if motif(1) >= ADAPTER_SUGGEST_RATE {
            suggestions.push(Suggestion {
                flags: vec!["--adapter-preset".into(), "nextera".into()],
                reason: format!("Nextera adapter motif in {:.2}% of reads", motif(1) * 100.0),
            });
        }
    }
    for (m, d) in detections.iter().enumerate() {
        if let (None, Some(cand)) = (&d.adapter, &d.candidate) {
            notes.push(format!(
                "R{}: over-represented sequence {cand} (not a known adapter); pass it with --adapter-r{} if it is one",
                m + 1,
                m + 1
            ));
        }
        if d.adapter.is_none() && c.samples[m].len() < detect::MIN_READS {
            notes.push(format!(
                "R{}: too few reads for adapter detection (need {})",
                m + 1,
                group(detect::MIN_READS as u64)
            ));
        }
    }

    // quality
    if let Some((m, cycle)) = low_tail_cycle {
        suggestions.push(Suggestion {
            flags: vec!["--cut-tail".into()],
            reason: format!("R{} mean quality drops below Q20 at cycle {cycle}", m + 1),
        });
    }
    if low_front {
        suggestions.push(Suggestion {
            flags: vec!["--cut-front".into()],
            reason: "mean quality of the first cycles is below Q20".into(),
        });
    }
    let median_len = hist_quantile(&st.mates[0].len_hist, 0.5).unwrap_or(0);
    if median_len > 1000 {
        warnings.push("reads are long (median > 1000 bp); this looks like long-read data, for which 'simpletrim long' is meant".into());
    }
    if st.mates[..mates]
        .iter()
        .any(|ms| hist_min_max(&ms.len_hist).is_some_and(|(a, b)| a != b))
    {
        notes.push("read lengths vary; the reads may have been trimmed already".into());
    }

    // suggested command
    let mut cmd = vec![
        "simpletrim".to_string(),
        "short".into(),
        "--r1".into(),
        args.r1.display().to_string(),
    ];
    if let Some(r2) = &args.r2 {
        cmd.extend(["--r2".into(), r2.display().to_string()]);
    }
    cmd.extend(["--o1".into(), trimmed_name(&args.r1)]);
    if let Some(r2) = &args.r2 {
        cmd.extend(["--o2".into(), trimmed_name(r2)]);
    }
    for s in &suggestions {
        cmd.extend(s.flags.iter().cloned());
    }
    let command = cmd
        .iter()
        .map(|a| shell_quote(a))
        .collect::<Vec<_>>()
        .join(" ");

    let insert_json = paired.then(|| {
        json!({
            "pairs_overlapping": ps.overlapped,
            "overlapping_fraction": round(if ps.pairs > 0 { ps.overlapped as f64 / ps.pairs as f64 } else { 0.0 }, 5),
            "adapter_read_through_fraction": round(if ps.pairs > 0 { ps.read_through as f64 / ps.pairs as f64 } else { 0.0 }, 5),
            "insert_size": {
                "median": hist_quantile(&ps.insert_hist, 0.5),
                "mode": hist_mode(&ps.insert_hist),
                "p10": hist_quantile(&ps.insert_hist, 0.1),
                "p90": hist_quantile(&ps.insert_hist, 0.9),
                "capped_at": MAX_INSERT,
            },
        })
    });
    let report = json!({
        "version": env!("CARGO_PKG_VERSION"),
        "status": "ok",
        "paired": paired,
        "inputs": { "r1": file_json(&args.r1), "r2": args.r2.as_deref().map(file_json) },
        "platform": header_json(&h1),
        "sample": {
            "records": total,
            "whole_file": c.whole_file,
            "estimated_total_records": c.estimated_total,
        },
        "pairing": paired.then(|| json!({
            "name_mismatches": ps.name_mismatches,
            "first_mismatch": ps.first_mismatch.as_ref().map(|(n, a, b)| json!({"pair": n, "r1": a, "r2": b})),
            "read_numbers": [h1.read_number, h2.as_ref().and_then(|h| h.read_number)],
            "count_mismatch": c.count_mismatch,
        })),
        "reads": if paired { json!({"r1": mate_json[0], "r2": mate_json[1]}) } else { json!({"r1": mate_json[0]}) },
        "overlap": insert_json,
        "warnings": warnings,
        "notes": notes,
        "suggestions": suggestions.iter().map(|s| json!({"flags": s.flags, "reason": s.reason})).collect::<Vec<_>>(),
        "suggested_command": command,
    });

    let text = render_text(args, &report, mates);
    emit(log, &text, &report)
}

fn pct(v: &Value) -> String {
    v.as_f64()
        .map_or("-".into(), |f| format!("{:.2}%", f * 100.0))
}

fn num(v: &Value) -> String {
    match v {
        Value::Number(n) if n.is_u64() => group(n.as_u64().unwrap_or(0)),
        Value::Number(n) => format!("{:.2}", n.as_f64().unwrap_or(0.0)),
        Value::Null => "-".into(),
        Value::String(s) => s.clone(),
        other => other.to_string(),
    }
}

fn render_text(args: &AnalyzeArgs, r: &Value, mates: usize) -> String {
    let mut t = String::new();
    let mut line = |s: String| {
        t.push_str(&s);
        t.push('\n');
    };
    let unit = if mates == 2 { "read pairs" } else { "reads" };
    let sample = &r["sample"];
    let scope = if sample["whole_file"].as_bool() == Some(true) {
        "whole file".to_string()
    } else {
        match sample["estimated_total_records"].as_u64() {
            Some(est) => format!("first part of ~{} {unit}", group(est)),
            None => "first part of the file".into(),
        }
    };
    line(format!(
        "simpletrim analyze: {} {unit} inspected ({scope})",
        num(&sample["records"])
    ));
    line(String::new());
    line("Input".into());
    for (m, path) in [Some(&args.r1), args.r2.as_ref()].into_iter().enumerate() {
        let Some(path) = path else { continue };
        let f = &r["inputs"][format!("r{}", m + 1)];
        let size = f["file_size"]
            .as_u64()
            .map(|b| format!(", {:.2} GB", b as f64 / 1e9))
            .unwrap_or_default();
        line(format!(
            "  R{}: {} ({}{size})",
            m + 1,
            path.display(),
            if f["compressed"].as_bool() == Some(true) {
                "gzip"
            } else {
                "plain"
            }
        ));
    }
    let p = &r["platform"];
    let mut platform = p["name"].as_str().unwrap_or("unknown").to_string();
    if let Some(model) = p["instrument_model"].as_str()
        && !platform.contains(model)
    {
        platform.push_str(&format!(" {model}"));
    }
    if let Some(id) = p["instrument_id"].as_str() {
        platform.push_str(&format!(" (instrument {id})"));
    }
    line(format!(
        "  Platform: {platform}; header format {}",
        num(&p["header_format"])
    ));
    line(format!(
        "  Two-colour chemistry: {}",
        if p["two_color"].as_bool() == Some(true) {
            "yes (poly-G tails expected)"
        } else {
            "no"
        }
    ));
    if mates == 2 {
        let pairing = &r["pairing"];
        let mism = pairing["name_mismatches"].as_u64().unwrap_or(0);
        line(format!(
            "  Pairing: {}; header read numbers R1={}, R2={}",
            if mism == 0 {
                "read names match in all pairs".to_string()
            } else {
                format!("{} pairs with different names", group(mism))
            },
            num(&pairing["read_numbers"][0]),
            num(&pairing["read_numbers"][1])
        ));
    }
    line(String::new());

    let reads = &r["reads"];
    let mate = |m: usize| &reads[format!("r{}", m + 1)];
    let row = |label: &str, f: &dyn Fn(&Value) -> String| {
        let mut s = format!("  {label:<30}");
        for m in 0..mates {
            s.push_str(&format!("{:>22}", f(mate(m))));
        }
        s
    };
    line(if mates == 2 {
        format!("Reads{:>47}{:>22}", "R1", "R2")
    } else {
        format!("Reads{:>47}", "R1")
    });
    line(row("reads", &|v| num(&v["reads"])));
    line(row("length min / median / max", &|v| {
        format!(
            "{} / {} / {}",
            num(&v["length"]["min"]),
            num(&v["length"]["median"]),
            num(&v["length"]["max"])
        )
    }));
    line(row("mean length", &|v| num(&v["length"]["mean"])));
    line(row("mean quality", &|v| num(&v["quality"]["mean"])));
    line(row("Q20 / Q30 bases", &|v| {
        format!(
            "{} / {}",
            pct(&v["quality"]["q20_rate"]),
            pct(&v["quality"]["q30_rate"])
        )
    }));
    line(row("quality encoding", &|v| num(&v["quality"]["encoding"])));
    line(row("distinct quality values", &|v| {
        format!(
            "{}{}",
            num(&v["quality"]["distinct_values"]),
            if v["quality"]["binned"].as_bool() == Some(true) {
                " (binned)"
            } else {
                ""
            }
        )
    }));
    line(row("first cycle < Q30 / < Q20", &|v| {
        format!(
            "{} / {}",
            num(&v["quality"]["first_cycle_below_q30"]),
            num(&v["quality"]["first_cycle_below_q20"])
        )
    }));
    line(row("GC content", &|v| pct(&v["gc_content"])));
    line(row("N bases", &|v| pct(&v["n_rate"])));
    line(row("poly-G tails (>= 10 bp)", &|v| {
        pct(&v["poly_g"]["fraction"])
    }));
    for (key, name, _) in MOTIFS {
        line(row(name, &|v| pct(&v["motifs"][key]["fraction"])));
    }
    line(row("composition bias cycles", &|v| {
        let cycles: Vec<String> = v["composition_bias_cycles"]
            .as_array()
            .map(|a| a.iter().map(|c| c.to_string()).collect())
            .unwrap_or_default();
        if cycles.is_empty() {
            "none".into()
        } else {
            cycles.join(",")
        }
    }));
    line(String::new());

    if let Some(ov) = r["overlap"].as_object() {
        line("Read overlap (paired-end)".into());
        line(format!(
            "  pairs overlapping: {}; with adapter read-through: {}",
            pct(&ov["overlapping_fraction"]),
            pct(&ov["adapter_read_through_fraction"])
        ));
        let is = &ov["insert_size"];
        line(format!(
            "  insert size of overlapping pairs: median {}, mode {}, 10-90%: {}-{}",
            num(&is["median"]),
            num(&is["mode"]),
            num(&is["p10"]),
            num(&is["p90"])
        ));
        line(String::new());
    }

    line("Adapters".into());
    for m in 0..mates {
        let a = &mate(m)["adapter"];
        let desc = match a["detected"].as_object() {
            Some(d) => format!(
                "{}{} - in {} of {} sampled reads",
                d["sequence"].as_str().unwrap_or(""),
                d["name"]
                    .as_str()
                    .map(|n| format!(" ({n})"))
                    .unwrap_or_default(),
                pct(&a["reads_with_adapter_fraction"]),
                num(&a["sample_reads"])
            ),
            None => "none detected".into(),
        };
        line(format!("  R{}: {desc}", m + 1));
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
    fn histogram_helpers() {
        let h = [0, 2, 0, 5, 1];
        assert_eq!(hist_quantile(&h, 0.5), Some(3));
        assert_eq!(hist_quantile(&h, 0.1), Some(1));
        assert_eq!(hist_mode(&h), Some(3));
        assert_eq!(hist_min_max(&h), Some((1, 4)));
        assert_eq!(hist_quantile(&[0, 0], 0.5), None);
    }

    #[test]
    fn output_names() {
        assert_eq!(
            trimmed_name(Path::new("dir/s_R1.fastq.gz")),
            "dir/s_R1.trimmed.fastq.gz"
        );
        assert_eq!(trimmed_name(Path::new("s.FQ")), "s.trimmed.FQ");
        assert_eq!(trimmed_name(Path::new("reads")), "reads.trimmed.fastq.gz");
    }

    #[test]
    fn quoting() {
        assert_eq!(shell_quote("a/b.fq.gz"), "a/b.fq.gz");
        assert_eq!(shell_quote("my file.fq"), "\"my file.fq\"");
    }

    #[test]
    fn encodings() {
        let mut h = vec![0u64; 256];
        h[b'#' as usize] = 1;
        h[b'F' as usize] = 1;
        assert_eq!(quality_encoding(&h).0, "phred+33");
        let mut h = vec![0u64; 256];
        h[b'B' as usize] = 1;
        h[b'h' as usize] = 1;
        assert_eq!(quality_encoding(&h).0, "phred+64 (likely)");
    }
}
