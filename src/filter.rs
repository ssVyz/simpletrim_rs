//! Per-read quality and length filters (port of fastp's `Filter::passFilter`).

use serde::Serialize;

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
}
