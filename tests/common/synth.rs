//! Synthetic long reads with a known structure, for tests and calibration.
//!
//! Reads look like nanopore reads: an adapter structure at the start (often missing its first
//! few bases), one insert (two for chimeras, joined by end and start adapter), and sometimes
//! a truncated adapter at the end. Errors are mostly insertions and deletions, and erroneous
//! bases get low qualities. Adapters are random sequences, so nothing relies on real kit
//! sequences.
//!
//! This file is used by the unit tests (included with `#[path]`) and the integration tests.
#![allow(dead_code)]

/// Small deterministic random number generator (SplitMix64).
pub struct Rng(u64);

impl Rng {
    pub fn new(seed: u64) -> Self {
        Rng(seed)
    }

    pub fn next_u64(&mut self) -> u64 {
        self.0 = self.0.wrapping_add(0x9E37_79B9_7F4A_7C15);
        let mut z = self.0;
        z = (z ^ (z >> 30)).wrapping_mul(0xBF58_476D_1CE4_E5B9);
        z = (z ^ (z >> 27)).wrapping_mul(0x94D0_49BB_1331_11EB);
        z ^ (z >> 31)
    }

    pub fn below(&mut self, n: usize) -> usize {
        (self.next_u64() % n as u64) as usize
    }

    /// Uniform integer in `lo..=hi`.
    pub fn range(&mut self, lo: usize, hi: usize) -> usize {
        lo + self.below(hi - lo + 1)
    }

    pub fn uniform(&mut self) -> f64 {
        (self.next_u64() >> 11) as f64 / (1u64 << 53) as f64
    }

    pub fn chance(&mut self, p: f64) -> bool {
        self.uniform() < p
    }

    pub fn seq(&mut self, len: usize) -> Vec<u8> {
        (0..len).map(|_| b"ACGT"[self.below(4)]).collect()
    }
}

/// Per-base error rates.
#[derive(Clone, Copy, Debug)]
pub struct Errors {
    pub sub: f64,
    pub ins: f64,
    pub del: f64,
}

impl Errors {
    pub const NONE: Errors = Errors {
        sub: 0.0,
        ins: 0.0,
        del: 0.0,
    };

    /// A total error rate split into 40% deletions, 30% insertions and 30% substitutions.
    pub fn nanopore(total: f64) -> Errors {
        Errors {
            sub: 0.3 * total,
            ins: 0.3 * total,
            del: 0.4 * total,
        }
    }
}

/// Appends a copy of `src` with random errors. Correct bases get a quality in `good`,
/// erroneous ones a quality between 2 and 7.
pub fn mutate_into(
    rng: &mut Rng,
    src: &[u8],
    e: Errors,
    good: (u8, u8),
    seq: &mut Vec<u8>,
    qual: &mut Vec<u8>,
) {
    let bad_q = |rng: &mut Rng| 33 + rng.range(2, 7) as u8;
    for &b in src {
        if rng.chance(e.ins) {
            seq.push(b"ACGT"[rng.below(4)]);
            qual.push(bad_q(rng));
        }
        if rng.chance(e.del) {
            continue;
        }
        if rng.chance(e.sub) {
            let other = b"ACGT"
                .iter()
                .copied()
                .filter(|&x| x != b)
                .nth(rng.below(3))
                .unwrap();
            seq.push(other);
            qual.push(bad_q(rng));
        } else {
            seq.push(b);
            qual.push(33 + rng.range(good.0 as usize, good.1 as usize) as u8);
        }
    }
}

/// Adapter structures of a synthetic library.
pub struct Library {
    pub start_adapter: Vec<u8>,
    pub end_adapter: Vec<u8>,
}

impl Library {
    /// Random adapters: a 70 bp start structure (adapter plus barcode) and a 40 bp end adapter.
    pub fn random(rng: &mut Rng) -> Self {
        Library {
            start_adapter: rng.seq(70),
            end_adapter: rng.seq(40),
        }
    }
}

#[derive(Clone, Debug)]
pub struct Params {
    /// Insert lengths are log-uniform in this range.
    pub insert_len: (usize, usize),
    pub errors: Errors,
    /// Quality range of correct bases.
    pub good_quality: (u8, u8),
    pub start_adapter_rate: f64,
    pub end_adapter_rate: f64,
    pub chimera_rate: f64,
    /// Fraction of reads (with inserts of at least 1500 bp) that get a 300 bp stretch of very
    /// poor sequence inside the insert.
    pub low_quality_rate: f64,
}

impl Default for Params {
    fn default() -> Self {
        Params {
            insert_len: (300, 5000),
            errors: Errors::nanopore(0.06),
            good_quality: (12, 30),
            start_adapter_rate: 0.9,
            end_adapter_rate: 0.5,
            chimera_rate: 0.0,
            low_quality_rate: 0.0,
        }
    }
}

/// A simulated read and where its parts are.
#[derive(Clone, Debug)]
pub struct SimRead {
    pub name: String,
    pub seq: Vec<u8>,
    pub qual: Vec<u8>,
    /// `[start, end)` of each insert in the read.
    pub inserts: Vec<(usize, usize)>,
    /// `[start, end)` of poor stretches inside inserts.
    pub low_quality: Vec<(usize, usize)>,
    pub has_start_adapter: bool,
    pub has_end_adapter: bool,
}

pub fn simulate(rng: &mut Rng, lib: &Library, p: &Params, index: usize) -> SimRead {
    let mut seq = Vec::new();
    let mut qual = Vec::new();
    let e = p.errors;
    let good = p.good_quality;
    let has_start_adapter = rng.chance(p.start_adapter_rate);
    if has_start_adapter {
        // the first few adapter bases are often lost
        let missing = rng.below(8);
        mutate_into(
            rng,
            &lib.start_adapter[missing..],
            e,
            good,
            &mut seq,
            &mut qual,
        );
    }
    let parts = if rng.chance(p.chimera_rate) { 2 } else { 1 };
    let mut inserts = Vec::new();
    let mut low_quality = Vec::new();
    for k in 0..parts {
        if k > 0 {
            mutate_into(rng, &lib.end_adapter, e, good, &mut seq, &mut qual);
            mutate_into(rng, &lib.start_adapter, e, good, &mut seq, &mut qual);
        }
        let (lo, hi) = (p.insert_len.0 as f64, p.insert_len.1 as f64);
        let len = (lo * (hi / lo).powf(rng.uniform())).round() as usize;
        let insert = rng.seq(len);
        let start = seq.len();
        if len >= 1500 && rng.chance(p.low_quality_rate) {
            let bad_start = rng.range(500, len - 800);
            mutate_into(rng, &insert[..bad_start], e, good, &mut seq, &mut qual);
            let s = seq.len();
            mutate_into(
                rng,
                &insert[bad_start..bad_start + 300],
                Errors::nanopore(0.35),
                (2, 6),
                &mut seq,
                &mut qual,
            );
            low_quality.push((s, seq.len()));
            mutate_into(
                rng,
                &insert[bad_start + 300..],
                e,
                good,
                &mut seq,
                &mut qual,
            );
        } else {
            mutate_into(rng, &insert, e, good, &mut seq, &mut qual);
        }
        inserts.push((start, seq.len()));
    }
    let has_end_adapter = rng.chance(p.end_adapter_rate);
    if has_end_adapter {
        // the read often ends inside the adapter
        let keep = rng.range(12, lib.end_adapter.len());
        mutate_into(rng, &lib.end_adapter[..keep], e, good, &mut seq, &mut qual);
    }
    SimRead {
        name: format!("read{index}"),
        seq,
        qual,
        inserts,
        low_quality,
        has_start_adapter,
        has_end_adapter,
    }
}

/// Simulates `n` reads.
pub fn simulate_many(seed: u64, lib: &Library, p: &Params, n: usize) -> Vec<SimRead> {
    let mut rng = Rng::new(seed);
    (0..n).map(|i| simulate(&mut rng, lib, p, i)).collect()
}

/// FASTQ text of the reads (names as given, no comments).
pub fn to_fastq(reads: &[SimRead]) -> String {
    let mut out = String::new();
    for r in reads {
        out.push('@');
        out.push_str(&r.name);
        out.push('\n');
        out.push_str(std::str::from_utf8(&r.seq).unwrap());
        out.push_str("\n+\n");
        out.push_str(std::str::from_utf8(&r.qual).unwrap());
        out.push('\n');
    }
    out
}
