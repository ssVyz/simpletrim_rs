//! Per-read quality and length filters (port of fastp's and fastplong's `Filter::passFilter`).

use serde::Serialize;

use crate::quality::mean_error_quality;

/// fastp's default quality filter: discard reads with more than 40% of bases below Q15,
/// or with more than 5 `N` bases. An average quality requirement is off by default.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct QualityFilter {
    pub qualified_quality: u8,
    pub unqualified_percent_limit: u32,
    pub n_base_limit: usize,
    pub average_quality: u32,
}

/// Reads shorter than `min_len` or longer than `max_len` (if non-zero) are discarded.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct LengthFilter {
    pub min_len: usize,
    pub max_len: usize,
}

/// What happened to a read.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Fate {
    Pass,
    TooShort,
    TooLong,
    LowQuality,
    TooManyN,
    AdapterDimer,
}

/// Applies the filters to a (trimmed) read. Empty reads always fail as too short.
pub fn evaluate(
    seq: &[u8],
    qual: &[u8],
    quality: Option<&QualityFilter>,
    length: Option<&LengthFilter>,
) -> Fate {
    let rlen = seq.len();
    if rlen == 0 {
        return Fate::TooShort;
    }
    if let Some(qf) = quality {
        let threshold = 33 + qf.qualified_quality;
        let mut low_qual = 0usize;
        let mut n_bases = 0usize;
        let mut total_qual = 0i64;
        for (&b, &q) in seq.iter().zip(qual) {
            total_qual += q as i64 - 33;
            low_qual += (q < threshold) as usize;
            n_bases += (b == b'N') as usize;
        }
        if low_qual as f64 > (qf.unqualified_percent_limit as usize * rlen) as f64 / 100.0 {
            return Fate::LowQuality;
        }
        if qf.average_quality > 0 && total_qual / (rlen as i64) < qf.average_quality as i64 {
            return Fate::LowQuality;
        }
        if n_bases > qf.n_base_limit {
            return Fate::TooManyN;
        }
    }
    if let Some(lf) = length {
        if rlen < lf.min_len {
            return Fate::TooShort;
        }
        if lf.max_len > 0 && rlen > lf.max_len {
            return Fate::TooLong;
        }
    }
    Fate::Pass
}

/// fastplong's quality filter: discard reads with more than 40% of bases below Q15 or more than
/// 10% `N` bases. Optional: a maximum number of `N` bases and a minimum mean quality, which is
/// computed from error probabilities here (fastplong averages phred values).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongQualityFilter {
    pub qualified_quality: u8,
    pub unqualified_percent_limit: u32,
    pub n_percent_limit: u32,
    pub n_base_limit: Option<usize>,
    pub average_quality: u8,
}

/// All filters for long reads; `None` disables a filter.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct LongFilter {
    pub quality: Option<LongQualityFilter>,
    pub length: Option<LengthFilter>,
    /// Allowed GC content in percent.
    pub gc: Option<(f64, f64)>,
    /// Minimum fraction of bases that differ from the next base.
    pub complexity: Option<f64>,
}

/// What happened to a long read (or a piece of one).
#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum LongFate {
    Pass,
    TooShort,
    TooLong,
    LowQuality,
    TooManyN,
    GcContent,
    LowComplexity,
    Chimeric,
}

impl LongFate {
    /// Tag appended to the read name in the failed-reads output (fastplong's names).
    pub fn tag(self) -> &'static str {
        match self {
            LongFate::Pass => "passed",
            LongFate::TooShort => "failed_too_short",
            LongFate::TooLong => "failed_too_long",
            LongFate::LowQuality => "failed_quality_filter",
            LongFate::TooManyN => "failed_too_many_n_bases",
            LongFate::GcContent => "failed_gc_content_filter",
            LongFate::LowComplexity => "failed_low_complexity",
            LongFate::Chimeric => "failed_chimeric_read",
        }
    }
}

/// Applies the long-read filters in fastplong's order: quality, length, GC content,
/// complexity. Empty reads always fail as too short.
pub fn evaluate_long(seq: &[u8], qual: &[u8], f: &LongFilter) -> LongFate {
    let rlen = seq.len();
    if rlen == 0 {
        return LongFate::TooShort;
    }
    if let Some(qf) = &f.quality {
        let threshold = 33 + qf.qualified_quality;
        let low_qual = qual.iter().filter(|&&q| q < threshold).count();
        let n_bases = seq.iter().filter(|&&b| b == b'N').count();
        if low_qual as f64 > (qf.unqualified_percent_limit as usize * rlen) as f64 / 100.0 {
            return LongFate::LowQuality;
        }
        if qf.average_quality > 0
            && mean_error_quality(qual).is_some_and(|q| q < qf.average_quality as f64)
        {
            return LongFate::LowQuality;
        }
        if n_bases * 100 > rlen * qf.n_percent_limit as usize {
            return LongFate::TooManyN;
        }
        if qf.n_base_limit.is_some_and(|limit| n_bases > limit) {
            return LongFate::TooManyN;
        }
    }
    if let Some(lf) = &f.length {
        if rlen < lf.min_len {
            return LongFate::TooShort;
        }
        if lf.max_len > 0 && rlen > lf.max_len {
            return LongFate::TooLong;
        }
    }
    if let Some((min, max)) = f.gc {
        let gc = seq
            .iter()
            .filter(|&&b| matches!(b, b'G' | b'C' | b'g' | b'c'))
            .count();
        let percent = gc as f64 * 100.0 / rlen as f64;
        if percent < min || percent > max {
            return LongFate::GcContent;
        }
    }
    if let Some(threshold) = f.complexity {
        let diff = seq.windows(2).filter(|w| w[0] != w[1]).count();
        if rlen <= 1 || (diff as f64) / ((rlen - 1) as f64) < threshold {
            return LongFate::LowComplexity;
        }
    }
    LongFate::Pass
}

#[cfg(test)]
mod tests {
    use super::*;

    const QF: QualityFilter = QualityFilter {
        qualified_quality: 15,
        unqualified_percent_limit: 40,
        n_base_limit: 5,
        average_quality: 0,
    };
    const LF: LengthFilter = LengthFilter {
        min_len: 15,
        max_len: 0,
    };

    #[test]
    fn passes_good_read() {
        let seq = b"ACGTACGTACGTACGTACGT";
        let qual = b"IIIIIIIIIIIIIIIIIIII";
        assert_eq!(evaluate(seq, qual, Some(&QF), Some(&LF)), Fate::Pass);
    }

    #[test]
    fn low_quality_and_n() {
        let seq = b"ACGTACGTACGTACGTACGT";
        assert_eq!(
            evaluate(seq, b"########IIIIIIIIIIII", Some(&QF), Some(&LF)),
            Fate::Pass
        );
        assert_eq!(
            evaluate(seq, b"#########IIIIIIIIIII", Some(&QF), Some(&LF)),
            Fate::LowQuality
        );
        assert_eq!(
            evaluate(
                b"ACGTNNNNNNACGTACGTAC",
                b"IIIIIIIIIIIIIIIIIIII",
                Some(&QF),
                Some(&LF)
            ),
            Fate::TooManyN
        );
    }

    #[test]
    fn length_limits() {
        let lf = LengthFilter {
            min_len: 5,
            max_len: 8,
        };
        assert_eq!(evaluate(b"ACGT", b"IIII", None, Some(&lf)), Fate::TooShort);
        assert_eq!(
            evaluate(b"ACGTACGTA", b"IIIIIIIII", None, Some(&lf)),
            Fate::TooLong
        );
        assert_eq!(evaluate(b"", b"", None, None), Fate::TooShort);
        assert_eq!(evaluate(b"ACGTAC", b"IIIIII", None, Some(&lf)), Fate::Pass);
    }

    #[test]
    fn average_quality() {
        let qf = QualityFilter {
            average_quality: 30,
            ..QF
        };
        // Q20 everywhere: all bases qualified, but the mean is below 30
        assert_eq!(
            evaluate(b"ACGTACGTAC", b"5555555555", Some(&qf), None),
            Fate::LowQuality
        );
    }

    const LONG: LongFilter = LongFilter {
        quality: Some(LongQualityFilter {
            qualified_quality: 15,
            unqualified_percent_limit: 40,
            n_percent_limit: 10,
            n_base_limit: None,
            average_quality: 0,
        }),
        length: Some(LengthFilter {
            min_len: 5,
            max_len: 0,
        }),
        gc: None,
        complexity: None,
    };

    #[test]
    fn long_read_filters() {
        let seq = b"ACGTACGTAC";
        assert_eq!(evaluate_long(seq, b"IIIIIIIIII", &LONG), LongFate::Pass);
        // N percent: 1 of 10 is allowed, 2 are not
        assert_eq!(
            evaluate_long(b"ACGTNCGTAC", b"IIIIIIIIII", &LONG),
            LongFate::Pass
        );
        assert_eq!(
            evaluate_long(b"NCGTNCGTAC", b"IIIIIIIIII", &LONG),
            LongFate::TooManyN
        );
        assert_eq!(evaluate_long(b"ACGT", b"IIII", &LONG), LongFate::TooShort);
        // mean quality from error probabilities: Q40 x 9 and Q2 x 1 is about Q12
        let mean = LongFilter {
            quality: Some(LongQualityFilter {
                average_quality: 13,
                ..LONG.quality.unwrap()
            }),
            ..LONG
        };
        assert_eq!(
            evaluate_long(seq, b"IIIII#IIII", &mean),
            LongFate::LowQuality
        );
        let gc = LongFilter {
            gc: Some((30.0, 60.0)),
            ..LONG
        };
        assert_eq!(
            evaluate_long(b"GCGCGCGCGC", b"IIIIIIIIII", &gc),
            LongFate::GcContent
        );
        assert_eq!(evaluate_long(seq, b"IIIIIIIIII", &gc), LongFate::Pass);
        let cx = LongFilter {
            complexity: Some(0.3),
            ..LONG
        };
        assert_eq!(
            evaluate_long(b"AAAAAAAAAC", b"IIIIIIIIII", &cx),
            LongFate::LowComplexity
        );
        assert_eq!(evaluate_long(seq, b"IIIIIIIIII", &cx), LongFate::Pass);
    }
}
