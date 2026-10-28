//! Counters collected while trimming, merged across chunks.

use std::collections::BTreeMap;

use serde::Serialize;
use serde_json::{Value, json};

use crate::filter::{Fate, LongFate};
use crate::pipeline::Merge;
use crate::quality::{error_sum, phred_of_error};

/// Adds every field of `other` to `self`.
macro_rules! impl_merge {
    ($ty:ty { $($field:ident),* $(,)? }) => {
        impl $ty {
            pub fn merge(&mut self, other: &$ty) {
                $( self.$field.merge_value(&other.$field); )*
            }
        }
    };
}

trait MergeValue {
    fn merge_value(&mut self, other: &Self);
}

impl MergeValue for u64 {
    fn merge_value(&mut self, other: &u64) {
        *self += other;
    }
}

/// Base-level statistics of a set of reads.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct BaseStats {
    pub reads: u64,
    pub bases: u64,
    pub q20_bases: u64,
    pub q30_bases: u64,
    pub gc_bases: u64,
    pub n_bases: u64,
}

impl_merge!(BaseStats {
    reads,
    bases,
    q20_bases,
    q30_bases,
    gc_bases,
    n_bases
});

impl BaseStats {
    pub fn add(&mut self, seq: &[u8], qual: &[u8]) {
        self.reads += 1;
        self.bases += seq.len() as u64;
        for (&b, &q) in seq.iter().zip(qual) {
            self.q20_bases += (q >= 33 + 20) as u64;
            self.q30_bases += (q >= 33 + 30) as u64;
            self.gc_bases += matches!(b, b'G' | b'C' | b'g' | b'c') as u64;
            self.n_bases += matches!(b, b'N' | b'n') as u64;
        }
    }

    fn rate(&self, n: u64) -> f64 {
        if self.bases == 0 {
            0.0
        } else {
            n as f64 / self.bases as f64
        }
    }

    pub fn mean_length(&self) -> f64 {
        if self.reads == 0 {
            0.0
        } else {
            self.bases as f64 / self.reads as f64
        }
    }

    pub fn q30_rate(&self) -> f64 {
        self.rate(self.q30_bases)
    }

    pub fn gc_content(&self) -> f64 {
        self.rate(self.gc_bases)
    }

    /// JSON form including derived rates.
    pub fn to_json(self) -> Value {
        json!({
            "reads": self.reads,
            "bases": self.bases,
            "mean_length": round(self.mean_length(), 2),
            "q20_bases": self.q20_bases,
            "q30_bases": self.q30_bases,
            "q20_rate": round(self.rate(self.q20_bases), 5),
            "q30_rate": round(self.q30_rate(), 5),
            "gc_content": round(self.gc_content(), 5),
            "n_bases": self.n_bases,
        })
    }
}

pub fn round(x: f64, digits: i32) -> f64 {
    let f = 10f64.powi(digits);
    (x * f).round() / f
}

/// Reads touched by a trimming step and bases removed by it.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct Trimmed {
    pub reads: u64,
    pub bases: u64,
}

impl_merge!(Trimmed { reads, bases });

impl MergeValue for Trimmed {
    fn merge_value(&mut self, other: &Trimmed) {
        self.merge(other);
    }
}

impl Trimmed {
    /// Counts a read if `bases > 0`.
    pub fn add(&mut self, bases: usize) {
        if bases > 0 {
            self.reads += 1;
            self.bases += bases as u64;
        }
    }
}

/// Bases removed per trimming step for one mate.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct TrimStats {
    pub fixed: Trimmed,
    pub quality_window: Trimmed,
    pub poly_g: Trimmed,
    pub adapter: Trimmed,
    pub poly_x: Trimmed,
}

impl_merge!(TrimStats {
    fixed,
    quality_window,
    poly_g,
    adapter,
    poly_x
});

/// How many reads of one mate ended with each [`Fate`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct FateCounts {
    pub passed: u64,
    pub too_short: u64,
    pub too_long: u64,
    pub low_quality: u64,
    pub too_many_n: u64,
    pub adapter_dimer: u64,
}

impl_merge!(FateCounts {
    passed,
    too_short,
    too_long,
    low_quality,
    too_many_n,
    adapter_dimer
});

impl FateCounts {
    pub fn add(&mut self, fate: Fate) {
        *match fate {
            Fate::Pass => &mut self.passed,
            Fate::TooShort => &mut self.too_short,
            Fate::TooLong => &mut self.too_long,
            Fate::LowQuality => &mut self.low_quality,
            Fate::TooManyN => &mut self.too_many_n,
            Fate::AdapterDimer => &mut self.adapter_dimer,
        } += 1;
    }
}

/// Outcomes of read pairs (paired-end only).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct PairCounts {
    pub pairs: u64,
    pub both_passed: u64,
    pub only_r1_passed: u64,
    pub only_r2_passed: u64,
    /// Includes adapter dimers.
    pub neither_passed: u64,
    pub adapter_dimers: u64,
    pub adapter_trimmed_by_overlap: u64,
    /// Pairs written to the paired outputs (both passed, or rescued in `keep-both` mode).
    pub pairs_written: u64,
    pub orphans_r1_written: u64,
    pub orphans_r2_written: u64,
    /// Failing mates written untrimmed in `keep-both` mode.
    pub rescued_r1_untrimmed: u64,
    pub rescued_r2_untrimmed: u64,
}

impl_merge!(PairCounts {
    pairs,
    both_passed,
    only_r1_passed,
    only_r2_passed,
    neither_passed,
    adapter_dimers,
    adapter_trimmed_by_overlap,
    pairs_written,
    orphans_r1_written,
    orphans_r2_written,
    rescued_r1_untrimmed,
    rescued_r2_untrimmed,
});

/// Everything counted during a `short` run. Index 0 is read 1, index 1 read 2.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct RunStats {
    pub input: [BaseStats; 2],
    pub output: [BaseStats; 2],
    pub orphans: BaseStats,
    pub trim: [TrimStats; 2],
    pub fates: [FateCounts; 2],
    pub pairs: PairCounts,
}

impl Merge for RunStats {
    fn merge(&mut self, other: &RunStats) {
        for m in 0..2 {
            self.input[m].merge(&other.input[m]);
            self.output[m].merge(&other.output[m]);
            self.trim[m].merge(&other.trim[m]);
            self.fates[m].merge(&other.fates[m]);
        }
        self.orphans.merge(&other.orphans);
        self.pairs.merge(&other.pairs);
    }
}

/// Distribution of read lengths (exact), for N50 and quantiles of long reads.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct LengthDist(BTreeMap<usize, u64>);

impl LengthDist {
    pub fn add(&mut self, len: usize) {
        *self.0.entry(len).or_insert(0) += 1;
    }

    pub fn merge(&mut self, other: &LengthDist) {
        for (&len, &n) in &other.0 {
            *self.0.entry(len).or_insert(0) += n;
        }
    }

    pub fn reads(&self) -> u64 {
        self.0.values().sum()
    }

    pub fn bases(&self) -> u64 {
        self.0.iter().map(|(&len, &n)| len as u64 * n).sum()
    }

    pub fn min(&self) -> Option<usize> {
        self.0.keys().next().copied()
    }

    pub fn max(&self) -> Option<usize> {
        self.0.keys().next_back().copied()
    }

    /// Length at quantile `q` of the reads.
    pub fn quantile(&self, q: f64) -> Option<usize> {
        let target = ((q * self.reads() as f64).ceil() as u64).max(1);
        let mut acc = 0;
        for (&len, &n) in &self.0 {
            acc += n;
            if acc >= target {
                return Some(len);
            }
        }
        None
    }

    /// The largest length L such that reads of at least L bases hold `fraction` of all bases
    /// (N50 for 0.5, N90 for 0.9).
    pub fn nx(&self, fraction: f64) -> Option<usize> {
        let target = fraction * self.bases() as f64;
        let mut acc = 0u64;
        for (&len, &n) in self.0.iter().rev() {
            acc += len as u64 * n;
            if acc as f64 >= target && acc > 0 {
                return Some(len);
            }
        }
        None
    }

    /// Bases in reads of at least `min_len` bases.
    pub fn bases_at_least(&self, min_len: usize) -> u64 {
        self.0
            .range(min_len..)
            .map(|(&len, &n)| len as u64 * n)
            .sum()
    }

    /// Lengths summary as JSON.
    pub fn to_json(&self) -> Value {
        let reads = self.reads();
        json!({
            "min": self.min(),
            "max": self.max(),
            "mean": round(if reads > 0 { self.bases() as f64 / reads as f64 } else { 0.0 }, 2),
            "median": self.quantile(0.5),
            "n50": self.nx(0.5),
            "n90": self.nx(0.9),
        })
    }
}

/// Highest read quality class; reads with a mean quality of at least this are counted in it.
const MAX_READ_QUALITY: usize = 60;
/// Read quality classes reported in JSON.
const READ_QUALITY_LEVELS: [usize; 5] = [7, 10, 12, 15, 20];

/// Statistics of a set of long reads. Quality means average error probabilities.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct ReadSetStats {
    pub base: BaseStats,
    pub lengths: LengthDist,
    /// Sum of the error probabilities of all bases.
    pub error_sum: f64,
    /// Number of reads per mean read quality, rounded down (index 60: 60 or more).
    pub read_quality: Vec<u64>,
}

impl ReadSetStats {
    pub fn add(&mut self, seq: &[u8], qual: &[u8]) {
        self.base.add(seq, qual);
        self.lengths.add(seq.len());
        if seq.is_empty() {
            return;
        }
        let errors = error_sum(qual);
        self.error_sum += errors;
        let q = phred_of_error(errors / seq.len() as f64);
        if self.read_quality.is_empty() {
            self.read_quality = vec![0; MAX_READ_QUALITY + 1];
        }
        self.read_quality[(q.max(0.0) as usize).min(MAX_READ_QUALITY)] += 1;
    }

    pub fn merge(&mut self, other: &ReadSetStats) {
        self.base.merge(&other.base);
        self.lengths.merge(&other.lengths);
        self.error_sum += other.error_sum;
        if self.read_quality.len() < other.read_quality.len() {
            self.read_quality.resize(other.read_quality.len(), 0);
        }
        for (a, b) in self.read_quality.iter_mut().zip(&other.read_quality) {
            *a += b;
        }
    }

    /// Mean quality over all bases.
    pub fn mean_quality(&self) -> Option<f64> {
        (self.base.bases > 0).then(|| phred_of_error(self.error_sum / self.base.bases as f64))
    }

    /// Median of the mean read qualities (rounded down).
    pub fn median_read_quality(&self) -> Option<usize> {
        let total: u64 = self.read_quality.iter().sum();
        let target = total.div_ceil(2).max(1);
        let mut acc = 0;
        for (q, &n) in self.read_quality.iter().enumerate() {
            acc += n;
            if acc >= target {
                return Some(q);
            }
        }
        None
    }

    /// Number of reads with a mean quality of at least `q`.
    pub fn reads_with_quality_at_least(&self, q: usize) -> u64 {
        self.read_quality.iter().skip(q).sum()
    }

    pub fn to_json(&self) -> Value {
        let mut v = self.base.to_json();
        v["length"] = self.lengths.to_json();
        v["mean_quality"] = json!(self.mean_quality().map(|q| round(q, 2)));
        v["median_read_quality"] = json!(self.median_read_quality());
        let levels: serde_json::Map<String, Value> = READ_QUALITY_LEVELS
            .iter()
            .map(|&q| (format!("q{q}"), json!(self.reads_with_quality_at_least(q))))
            .collect();
        v["reads_with_mean_quality_at_least"] = Value::Object(levels);
        v
    }
}

/// Bases removed per trimming step of `long`.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct LongTrimStats {
    pub fixed: Trimmed,
    pub quality_window: Trimmed,
    pub poly_x: Trimmed,
    pub adapter_start: Trimmed,
    pub adapter_end: Trimmed,
}

impl_merge!(LongTrimStats {
    fixed,
    quality_window,
    poly_x,
    adapter_start,
    adapter_end
});

/// Reads with adapters inside (chimeras).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct ChimeraCounts {
    pub reads: u64,
    pub reads_discarded: u64,
    /// Pieces kept after splitting (at least the minimum segment length).
    pub pieces: u64,
    /// Bases removed with the adapters and in pieces that were too short.
    pub bases_removed: u64,
}

impl_merge!(ChimeraCounts {
    reads,
    reads_discarded,
    pieces,
    bases_removed
});

/// Reads with low-quality stretches (`--low-quality-mode`).
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct LowQualityCounts {
    pub reads: u64,
    /// Pieces kept.
    pub pieces: u64,
    pub bases_removed: u64,
}

impl_merge!(LowQualityCounts {
    reads,
    pieces,
    bases_removed
});

/// How many records (reads or pieces of reads) ended with each [`LongFate`].
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq, Serialize)]
pub struct LongFateCounts {
    pub passed: u64,
    pub too_short: u64,
    pub too_long: u64,
    pub low_quality: u64,
    pub too_many_n: u64,
    pub gc_content: u64,
    pub low_complexity: u64,
    pub chimeric: u64,
}

impl_merge!(LongFateCounts {
    passed,
    too_short,
    too_long,
    low_quality,
    too_many_n,
    gc_content,
    low_complexity,
    chimeric
});

impl LongFateCounts {
    pub fn add(&mut self, fate: LongFate) {
        *match fate {
            LongFate::Pass => &mut self.passed,
            LongFate::TooShort => &mut self.too_short,
            LongFate::TooLong => &mut self.too_long,
            LongFate::LowQuality => &mut self.low_quality,
            LongFate::TooManyN => &mut self.too_many_n,
            LongFate::GcContent => &mut self.gc_content,
            LongFate::LowComplexity => &mut self.low_complexity,
            LongFate::Chimeric => &mut self.chimeric,
        } += 1;
    }
}

/// Everything counted during a `long` run.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct LongRunStats {
    pub input: ReadSetStats,
    pub output: ReadSetStats,
    pub trim: LongTrimStats,
    pub chimeras: ChimeraCounts,
    pub low_quality: LowQualityCounts,
    pub fates: LongFateCounts,
    /// Input reads with at least one record written.
    pub reads_with_output: u64,
    /// Input reads written as more than one record.
    pub reads_split: u64,
}

impl Merge for LongRunStats {
    fn merge(&mut self, other: &LongRunStats) {
        self.input.merge(&other.input);
        self.output.merge(&other.output);
        self.trim.merge(&other.trim);
        self.chimeras.merge(&other.chimeras);
        self.low_quality.merge(&other.low_quality);
        self.fates.merge(&other.fates);
        self.reads_with_output += other.reads_with_output;
        self.reads_split += other.reads_split;
    }
}

/// Formats an integer with thousands separators.
pub fn group(n: u64) -> String {
    let s = n.to_string();
    let mut out = String::with_capacity(s.len() + s.len() / 3);
    for (i, c) in s.chars().enumerate() {
        if i > 0 && (s.len() - i).is_multiple_of(3) {
            out.push(',');
        }
        out.push(c);
    }
    out
}

/// Percentage of `n` in `total`, as text.
pub fn percent(n: u64, total: u64) -> String {
    if total == 0 {
        "0.00%".into()
    } else {
        format!("{:.2}%", n as f64 * 100.0 / total as f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn base_stats() {
        let mut s = BaseStats::default();
        s.add(b"ACGN", b"5?IA");
        assert_eq!(s.reads, 1);
        assert_eq!(s.bases, 4);
        assert_eq!(s.q20_bases, 4);
        assert_eq!(s.q30_bases, 3);
        assert_eq!(s.gc_bases, 2);
        assert_eq!(s.n_bases, 1);
        let mut t = s;
        t.merge(&s);
        assert_eq!(t.bases, 8);
    }

    #[test]
    fn grouping() {
        assert_eq!(group(0), "0");
        assert_eq!(group(999), "999");
        assert_eq!(group(1000), "1,000");
        assert_eq!(group(1234567), "1,234,567");
        assert_eq!(percent(1, 4), "25.00%");
    }

    #[test]
    fn length_distribution() {
        let mut d = LengthDist::default();
        for len in [100, 200, 300, 400, 1000] {
            d.add(len);
        }
        assert_eq!(d.reads(), 5);
        assert_eq!(d.bases(), 2000);
        assert_eq!(d.quantile(0.5), Some(300));
        // 1000 alone is half of the bases
        assert_eq!(d.nx(0.5), Some(1000));
        assert_eq!(d.nx(0.9), Some(200));
        assert_eq!(d.bases_at_least(400), 1400);
        assert_eq!((d.min(), d.max()), (Some(100), Some(1000)));
        assert_eq!(LengthDist::default().nx(0.5), None);
    }

    #[test]
    fn read_set_quality() {
        let mut s = ReadSetStats::default();
        s.add(b"ACGT", b"++++"); // Q10
        s.add(b"ACGT", b"????"); // Q30
        assert!((s.mean_quality().unwrap() - 12.967).abs() < 0.01);
        assert_eq!(s.median_read_quality(), Some(10));
        assert_eq!(s.reads_with_quality_at_least(20), 1);
        let mut t = s.clone();
        t.merge(&s);
        assert_eq!(t.reads_with_quality_at_least(10), 4);
    }
}
