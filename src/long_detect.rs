//! Adapter auto-detection for long reads (port of fastplong's
//! `Evaluator::evalAdapterAndReadNum`).
//!
//! The 10-mers in the first 128 bases (read start) and in the last 128 bases (read end) of up
//! to 64K reads are counted. The most frequent informative 10-mer of each end is accepted as
//! an adapter seed if it is at least 100 times more frequent than the average 10-mer; it is
//! then extended base by base while one neighbouring 10-mer clearly dominates and sits at a
//! consistent position, up to 64 bases. Sequences longer than 16 bases are adapters.
//!
//! Difference to fastplong: the test that skips 10-mers with fewer than three base changes is
//! applied to the 10-mer; fastplong applies it to the 10-mer's count by mistake.

use crate::detect::{KEY_SPACE, KEYLEN, int2seq, seq2int};

/// fastplong reads up to 64K reads or 512 Mbp for detection.
pub const SAMPLE_READS: usize = 64 * 1024;
pub const SAMPLE_BASES: usize = 8192 * SAMPLE_READS;
/// Fewer reads than this are not enough to detect anything.
pub const MIN_READS: usize = 100;
/// Bases scanned at each read end.
const SCAN_LEN: usize = 128;
const FOLD_THRESHOLD: f64 = 100.0;
const MAX_LEN: usize = 64;
/// Detected sequences must be longer than this.
const MIN_LEN: usize = 16;
/// Reads checked for `U` (RNA) and `T` bases.
const RNA_CHECK_READS: usize = 100;

/// 10-mer counts at one read end.
struct EndCounts {
    counts: Vec<u32>,
    /// Sum of the positions (start: offset from the read start; end: distance from the read
    /// end) at which each 10-mer was seen.
    positions: Vec<u64>,
    total: u64,
}

impl EndCounts {
    fn new() -> Self {
        EndCounts {
            counts: vec![0; KEY_SPACE],
            positions: vec![0; KEY_SPACE],
            total: 0,
        }
    }

    fn add(&mut self, key: usize, position: usize) {
        self.counts[key] += 1;
        self.positions[key] += position as u64;
        self.total += 1;
    }
}

/// Result for one read end.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct EndDetection {
    pub adapter: Option<String>,
    /// A frequent sequence that is too short to be used (fastplong suggests passing it
    /// explicitly if it is an adapter).
    pub too_short: Option<String>,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Detection {
    pub start: EndDetection,
    pub end: EndDetection,
    pub reads: usize,
    /// The first reads contain `U` instead of `T` (direct RNA).
    pub rna: bool,
}

/// Collects 10-mer counts from reads; feed it with [`Detector::add`] while
/// [`Detector::wants_more`] is true.
pub struct Detector {
    start: EndCounts,
    end: EndCounts,
    shift_tail: usize,
    reads: usize,
    bases: usize,
    t_bases: u64,
    u_bases: u64,
}

impl Detector {
    /// `shift_tail` bases at the read end are ignored (fastplong: the fixed tail trimming,
    /// at least 1, because the last base is noisy).
    pub fn new(shift_tail: usize) -> Self {
        Detector {
            start: EndCounts::new(),
            end: EndCounts::new(),
            shift_tail: shift_tail.max(1),
            reads: 0,
            bases: 0,
            t_bases: 0,
            u_bases: 0,
        }
    }

    pub fn wants_more(&self) -> bool {
        self.reads < SAMPLE_READS && self.bases < SAMPLE_BASES
    }

    pub fn add(&mut self, seq: &[u8]) {
        if self.reads < RNA_CHECK_READS {
            self.t_bases += seq.iter().filter(|&&b| b == b'T').count() as u64;
            self.u_bases += seq.iter().filter(|&&b| b == b'U').count() as u64;
        }
        self.reads += 1;
        self.bases += seq.len();
        let len = seq.len();
        let Some(last) = len.checked_sub(KEYLEN + self.shift_tail) else {
            return;
        };
        let mut key = None;
        for pos in 0..=last.min(SCAN_LEN - 1) {
            key = seq2int(seq, pos, key);
            if let Some(k) = key {
                self.start.add(k, pos);
            }
        }
        let mut key = None;
        for pos in last.saturating_sub(SCAN_LEN)..=last {
            key = seq2int(seq, pos, key);
            if let Some(k) = key {
                self.end.add(k, len - pos);
            }
        }
    }

    pub fn finish(&self) -> Detection {
        let enough = self.reads >= MIN_READS;
        let rna = self.u_bases > 0 && self.t_bases == 0;
        let base = |b: u8| if rna && b == b'T' { b'U' } else { b };
        let end = |c: &EndCounts| {
            if !enough {
                return EndDetection::default();
            }
            let mut d = detect_end(c);
            for s in [&mut d.adapter, &mut d.too_short].into_iter().flatten() {
                *s = s.bytes().map(base).map(char::from).collect();
            }
            d
        };
        Detection {
            start: end(&self.start),
            end: end(&self.end),
            reads: self.reads,
            rna,
        }
    }
}

fn detect_end(c: &EndCounts) -> EndDetection {
    // number of distinct 10-mers, counted before AAAAAAAAAA is excluded (as fastplong)
    let distinct = c.counts.iter().filter(|&&n| n > 0).count() as f64;
    let Some(key) = top_key(&c.counts) else {
        return EndDetection::default();
    };
    let count = c.counts[key] as f64;
    if count <= 10.0 || count * distinct <= c.total as f64 * FOLD_THRESHOLD {
        return EndDetection::default();
    }
    let adapter = extend_key(key, c);
    if adapter.len() > MIN_LEN {
        EndDetection {
            adapter: Some(adapter),
            too_short: None,
        }
    } else {
        EndDetection {
            adapter: None,
            too_short: Some(adapter),
        }
    }
}

/// Base of a 10-mer at position `i` (0 = first base), in the A=0, T=1, C=2, G=3 encoding.
fn base_at(key: usize, i: usize) -> usize {
    (key >> ((KEYLEN - 1 - i) * 2)) & 3
}

/// The most frequent 10-mer that is not low-complexity, repetitive or GC-rich
/// (fastplong's `getTopKey`); the first one in key order on ties.
fn top_key(counts: &[u32]) -> Option<usize> {
    let mut top: Option<(usize, u32)> = None;
    for (k, &val) in counts.iter().enumerate() {
        if k == 0 || val == 0 || top.is_some_and(|(_, t)| val <= t) {
            continue;
        }
        let mut atcg = [0usize; 4];
        for i in 0..KEYLEN {
            atcg[base_at(k, i)] += 1;
        }
        let low_complexity =
            atcg.iter().any(|&n| n >= KEYLEN - 4) || atcg.iter().filter(|&&n| n == 0).count() >= 2;
        // first half equals second half
        let repetitive = k >> KEYLEN == k & ((1 << KEYLEN) - 1);
        let changes = (0..KEYLEN - 1)
            .filter(|&i| base_at(k, i) != base_at(k, i + 1))
            .count();
        let gc_rich = atcg[2] + atcg[3] >= KEYLEN - 2;
        let starts_with_gggg = k >> 12 == 0xff;
        if low_complexity || repetitive || changes < 3 || gc_rich || starts_with_gggg {
            continue;
        }
        top = Some((k, val));
    }
    top.map(|(k, _)| k)
}

/// Extends a seed 10-mer to the left, then to the right (fastplong's `extendKeyToAdapter`).
fn extend_key(key: usize, c: &EndCounts) -> String {
    const BASES: [u8; 4] = [b'A', b'T', b'C', b'G'];
    let mask = KEY_SPACE - 1;
    let mean_pos = |k: usize| c.positions[k] as f64 / c.counts[k] as f64;
    let mut adapter: std::collections::VecDeque<u8> = int2seq(key).into_iter().collect();
    for extending_left in [true, false] {
        let mut cur = key;
        while adapter.len() < MAX_LEN {
            let neighbour = |b: usize| {
                if extending_left {
                    (b << ((KEYLEN - 1) * 2)) | (cur >> 2)
                } else {
                    b | (mask & (cur << 2))
                }
            };
            let total: u64 = (0..4).map(|b| c.counts[neighbour(b)] as u64).sum();
            let next = (0..4).find(|&b| {
                let nk = neighbour(b);
                let n = c.counts[nk] as f64;
                if n == 0.0 {
                    return false;
                }
                let offset = mean_pos(nk) - mean_pos(cur);
                n / total as f64 >= 0.7
                    && n / c.counts[key] as f64 >= 0.5
                    && (-4.0..=2.0).contains(&offset)
            });
            let Some(b) = next else { break };
            cur = neighbour(b);
            if extending_left {
                adapter.push_front(BASES[b]);
            } else {
                adapter.push_back(BASES[b]);
            }
        }
    }
    adapter.into_iter().map(char::from).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::{Library, Params, Rng, simulate_many};

    fn detect(reads: &[Vec<u8>]) -> Detection {
        let mut d = Detector::new(1);
        for r in reads {
            d.add(r);
        }
        d.finish()
    }

    #[test]
    fn detects_synthetic_adapters() {
        let mut rng = Rng::new(1);
        let lib = Library::random(&mut rng);
        let reads: Vec<Vec<u8>> = simulate_many(2, &lib, &Params::default(), 3000)
            .into_iter()
            .map(|r| r.seq)
            .collect();
        let d = detect(&reads);
        let start = d.start.adapter.expect("start adapter");
        let end = d.end.adapter.expect("end adapter");
        let s = String::from_utf8(lib.start_adapter.clone()).unwrap();
        let e = String::from_utf8(lib.end_adapter.clone()).unwrap();
        // the detected sequences are (nearly exact) parts of the true adapters
        assert!(
            start.len() > 16 && s.contains(&start[2..start.len() - 2]),
            "{start} vs {s}"
        );
        assert!(
            end.len() > 16 && e.contains(&end[2..end.len() - 2]),
            "{end} vs {e}"
        );
        assert!(!d.rna);
    }

    #[test]
    fn nothing_in_reads_without_adapters() {
        let mut rng = Rng::new(3);
        let lib = Library::random(&mut rng);
        let p = Params {
            start_adapter_rate: 0.0,
            end_adapter_rate: 0.0,
            ..Params::default()
        };
        let reads: Vec<Vec<u8>> = simulate_many(4, &lib, &p, 3000)
            .into_iter()
            .map(|r| r.seq)
            .collect();
        let d = detect(&reads);
        assert_eq!(d.start, EndDetection::default());
        assert_eq!(d.end, EndDetection::default());
    }

    #[test]
    fn too_few_reads_and_rna() {
        let mut rng = Rng::new(5);
        let adapter = rng.seq(30);
        let make = |rng: &mut Rng| {
            let mut r = adapter.clone();
            r.extend(rng.seq(300));
            r
        };
        let few: Vec<Vec<u8>> = (0..50).map(|_| make(&mut rng)).collect();
        assert_eq!(detect(&few).start.adapter, None);
        let rna: Vec<Vec<u8>> = (0..500)
            .map(|_| {
                let mut r = make(&mut rng);
                r.iter_mut().filter(|b| **b == b'T').for_each(|b| *b = b'U');
                r
            })
            .collect();
        let d = detect(&rna);
        assert!(d.rna);
        let start = d.start.adapter.expect("adapter in RNA reads");
        assert!(start.contains('U') && !start.contains('T'));
    }

    #[test]
    fn top_key_skips_uninformative_kmers() {
        let mut counts = vec![0u32; KEY_SPACE];
        let code = |s: &[u8]| seq2int(s, 0, None).unwrap();
        counts[code(b"AAAAACCCCC")] = 1000; // low complexity
        counts[code(b"ACGTAACGTA")] = 900; // repetitive halves
        counts[code(b"AAAACCCGGG")] = 800; // only two base changes
        counts[code(b"GGGGCATCGC")] = 700; // GC-rich
        counts[code(b"ACTGATCAGT")] = 10;
        assert_eq!(top_key(&counts), Some(code(b"ACTGATCAGT")));
    }
}
