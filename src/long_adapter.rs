//! Adapter trimming for long reads: adapters at the read start and end (also when the read
//! starts or ends inside the adapter) and adapters inside reads (chimeras), all found by edit
//! distance (see [`crate::edit`]).
//!
//! Follows fastplong's `AdapterTrimmer`, with these differences:
//! - matches allow insertions and deletions everywhere (fastplong counts mismatches position by
//!   position and only checks the edit distance of one candidate position);
//! - all adapters of one read end are searched in the same window of the untrimmed read and
//!   the innermost match decides the cut (fastplong applies them one after another and moves
//!   its window after every cut);
//! - internal adapters use a stricter error rate and a minimum adapter length, because the
//!   whole read is searched and edit distance finds many more random matches than mismatch
//!   counting (see the calibration test in `edit.rs`).

use crate::edit::{Pattern, allowed, find_all, find_innermost};

/// Adapters shorter than this are not searched inside reads (too many random matches).
pub const MIN_INTERNAL_ADAPTER_LEN: usize = 20;

/// Matching parameters.
#[derive(Clone, Copy, Debug)]
pub struct Params {
    /// Allowed edits per aligned adapter base at the read ends (fastplong's
    /// `distance_threshold`).
    pub rate: f64,
    /// Minimum number of adapter bases for a match cut off by the read end.
    pub min_overlap: usize,
    /// Read bases searched at each end.
    pub window: usize,
    /// Extra bases removed next to every adapter match.
    pub extension: usize,
    /// Allowed edits per adapter base inside reads.
    pub internal_rate: f64,
}

/// Prepared adapter sequences.
pub struct LongAdapters {
    start: Vec<Pattern>,
    /// End adapters, reversed: the read end is searched backwards.
    end_reversed: Vec<Pattern>,
    internal: Vec<Pattern>,
    pub params: Params,
}

impl LongAdapters {
    /// `start` is searched at the read start, `end` at the read end, FASTA sequences at both
    /// ends; all of them (if long enough) inside reads.
    pub fn new(
        start: Option<&[u8]>,
        end: Option<&[u8]>,
        fasta: &[Vec<u8>],
        params: Params,
    ) -> Self {
        let mut at_start: Vec<&[u8]> = start.into_iter().collect();
        let mut at_end: Vec<&[u8]> = end.into_iter().collect();
        for f in fasta {
            at_start.push(f);
            at_end.push(f);
        }
        let mut internal: Vec<&[u8]> = Vec::new();
        for s in start
            .into_iter()
            .chain(end)
            .chain(fasta.iter().map(Vec::as_slice))
        {
            if s.len() >= MIN_INTERNAL_ADAPTER_LEN && !internal.contains(&s) {
                internal.push(s);
            }
        }
        LongAdapters {
            start: at_start.into_iter().map(Pattern::new).collect(),
            end_reversed: at_end
                .into_iter()
                .map(|s| Pattern::new(&s.iter().rev().copied().collect::<Vec<u8>>()))
                .collect(),
            internal: internal.into_iter().map(Pattern::new).collect(),
            params,
        }
    }

    /// Whether any adapter can be found inside reads.
    pub fn has_internal(&self) -> bool {
        !self.internal.is_empty()
    }

    /// Number of bases to remove from the start of `seq` (0 if no adapter is found).
    pub fn trim_start(&self, seq: &[u8]) -> usize {
        let window = &seq[..seq.len().min(self.params.window)];
        let end = self.innermost(&self.start, window);
        end.map_or(0, |e| (e + self.params.extension).min(seq.len()))
    }

    /// Number of bases to remove from the end of `seq` (0 if no adapter is found). `buf` is
    /// scratch space.
    pub fn trim_end(&self, seq: &[u8], buf: &mut Vec<u8>) -> usize {
        let w = seq.len().min(self.params.window);
        buf.clear();
        buf.extend(seq[seq.len() - w..].iter().rev());
        let end = self.innermost(&self.end_reversed, buf);
        end.map_or(0, |e| (e + self.params.extension).min(seq.len()))
    }

    /// Length of the innermost adapter match at the start of `text`.
    fn innermost(&self, patterns: &[Pattern], text: &[u8]) -> Option<usize> {
        patterns
            .iter()
            .filter_map(|p| {
                let min_overlap = self.params.min_overlap.min(p.len());
                find_innermost(p, text, self.params.rate, min_overlap).map(|h| h.end + 1)
            })
            .max()
    }

    /// Adapters inside `seq`: merged `[start, end)` ranges including the extension.
    pub fn internal_regions(&self, seq: &[u8]) -> Vec<(usize, usize)> {
        let ext = self.params.extension;
        let mut regions: Vec<(usize, usize)> = Vec::new();
        for p in &self.internal {
            let max_dist = allowed(self.params.internal_rate, p.len());
            for h in find_all(p, seq, max_dist) {
                let start = (h.end + 1).saturating_sub(p.len() + ext);
                let end = (h.end + 1 + ext).min(seq.len());
                regions.push((start, end));
            }
        }
        regions.sort_unstable();
        let mut merged: Vec<(usize, usize)> = Vec::with_capacity(regions.len());
        for (s, e) in regions {
            match merged.last_mut() {
                Some(last) if s <= last.1 => last.1 = last.1.max(e),
                _ => merged.push((s, e)),
            }
        }
        merged
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::synth::{Library, Params as SimParams, Rng, simulate_many};

    const PARAMS: Params = Params {
        rate: 0.25,
        min_overlap: 16,
        window: 200,
        extension: 10,
        internal_rate: 0.15,
    };

    #[test]
    fn exact_adapters_at_both_ends() {
        let start = b"ACTTCGTTCAGTTACGTATTGCTGGATCCAGTA";
        let end = b"GCAATACGTAACTGAACGAAGTACAGG";
        let ad = LongAdapters::new(Some(start), Some(end), &[], PARAMS);
        let insert = b"CATTTGACCAGATTACAGGATTTACCAGGGGTTTACCCATTACAGATACAGATGGACCAT".repeat(5);
        let mut read = start.to_vec();
        read.extend_from_slice(&insert);
        read.extend_from_slice(&end[..8]); // the read ends 8 bases into the adapter
        assert_eq!(ad.trim_start(&read), start.len() + 10);
        let mut buf = Vec::new();
        // too few adapter bases to be recognised
        assert_eq!(ad.trim_end(&read, &mut buf), 0);
        read.extend_from_slice(&end[8..20]);
        assert_eq!(ad.trim_end(&read, &mut buf), 20 + 10);
        // no adapter anywhere
        assert_eq!(ad.trim_start(&insert), 0);
        assert_eq!(ad.trim_end(&insert, &mut buf), 0);
    }

    #[test]
    fn internal_adapters_split_chimeras() {
        let start = b"ACTTCGTTCAGTTACGTATTGCTGGATCCAGTA";
        let ad = LongAdapters::new(Some(start), None, &[], PARAMS);
        let mut rng = Rng::new(5);
        let mut read = rng.seq(500);
        read.extend_from_slice(start);
        read.extend(rng.seq(400));
        let r = ad.internal_regions(&read);
        assert_eq!(r, vec![(490, 500 + start.len() + 10)]);
        // short adapters are not used inside reads
        let short = LongAdapters::new(Some(&start[..15]), None, &[], PARAMS);
        assert!(!short.has_internal());
    }

    /// Bases of adapter left at the read start / end, and insert bases removed, per read.
    fn trimming_errors(seed: u64, known: bool) -> (Vec<isize>, Vec<isize>) {
        let mut rng = Rng::new(seed);
        let lib = Library::random(&mut rng);
        let p = SimParams {
            insert_len: (300, 3000),
            ..SimParams::default()
        };
        let reads = simulate_many(seed, &lib, &p, 2000);
        let (start, end) = if known {
            (lib.start_adapter.clone(), lib.end_adapter.clone())
        } else {
            let mut det = crate::long_detect::Detector::new(1);
            for r in &reads {
                det.add(&r.seq);
            }
            let d = det.finish();
            (
                d.start
                    .adapter
                    .expect("start adapter detected")
                    .into_bytes(),
                d.end.adapter.expect("end adapter detected").into_bytes(),
            )
        };
        let ad = LongAdapters::new(Some(&start), Some(&end), &[], PARAMS);
        let mut buf = Vec::new();
        let (mut at_start, mut at_end) = (Vec::new(), Vec::new());
        for r in &reads {
            let (ins_start, ins_end) = r.inserts[0];
            let a = ad.trim_start(&r.seq);
            let b = r.seq.len() - ad.trim_end(&r.seq[a..], &mut buf);
            // positive: insert bases lost; negative: adapter bases left
            at_start.push(a as isize - ins_start as isize);
            at_end.push(ins_end as isize - b as isize);
        }
        (at_start, at_end)
    }

    fn summary(label: &str, errs: &[isize]) -> (f64, f64) {
        let n = errs.len() as f64;
        let left = errs.iter().filter(|&&e| e < -3).count() as f64 / n;
        let over = errs.iter().filter(|&&e| e > 20).count() as f64 / n;
        println!(
            "{label}: adapter left (> 3 bp) {:.2}%, insert lost (> 20 bp) {:.2}%",
            left * 100.0,
            over * 100.0
        );
        (left, over)
    }

    #[test]
    fn synthetic_reads_with_given_adapters() {
        let (s, e) = trimming_errors(11, true);
        let (left, over) = summary("start", &s);
        assert!(left < 0.02 && over < 0.01);
        let (left, over) = summary("end", &e);
        assert!(left < 0.02 && over < 0.01);
    }

    #[test]
    fn synthetic_reads_with_detected_adapters() {
        let (s, e) = trimming_errors(12, false);
        let (left, over) = summary("start", &s);
        assert!(left < 0.05 && over < 0.02);
        let (left, over) = summary("end", &e);
        assert!(left < 0.05 && over < 0.02);
    }
}
