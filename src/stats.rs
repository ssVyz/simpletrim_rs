//! Counters collected while trimming, merged across chunks.

use serde::Serialize;
use serde_json::{Value, json};

use crate::filter::Fate;

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
}

impl_merge!(TrimStats {
    fixed,
    quality_window,
    poly_g,
    adapter
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

impl RunStats {
    pub fn merge(&mut self, other: &RunStats) {
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
}
