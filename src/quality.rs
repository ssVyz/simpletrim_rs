//! Fixed front/tail trimming and sliding-window quality cutting (port of fastp's `trimAndCut`),
//! and detection of low-quality stretches inside reads (long reads).
//!
//! Window means are taken on a [`Scale`]: either the plain average of the phred values (fastp)
//! or the average of the base error probabilities, converted back to phred (the usual
//! definition for nanopore reads; it is dominated by the worst bases).

use std::sync::LazyLock;

/// One sliding window setting: window size and required mean phred quality.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Window {
    pub size: usize,
    pub quality: u8,
}

/// The three fastp sliding-window modes. `right` takes precedence over `tail`, as in fastp.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct QualityCut {
    pub front: Option<Window>,
    pub tail: Option<Window>,
    pub right: Option<Window>,
}

impl QualityCut {
    pub fn is_enabled(&self) -> bool {
        self.front.is_some() || self.tail.is_some() || self.right.is_some()
    }
}

/// How quality values are averaged over a window.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Averaging {
    /// Mean of the phred values (fastp, fastplong).
    Phred,
    /// Mean of the error probabilities `10^(-q/10)`, expressed as phred.
    ErrorProbability,
}

/// Per-base scores whose window sums are compared against a threshold: a window has a mean
/// quality of at least `Q` exactly if the sum of its scores is at least `size * threshold(Q)`.
/// Error probabilities are stored as negative fixed-point numbers, so all comparisons are exact
/// integer comparisons.
pub struct Scale {
    score: [i64; 256],
    averaging: Averaging,
}

/// Fixed-point unit of error probabilities.
const ERROR_UNIT: f64 = (1u64 << 30) as f64;

fn error_score(phred: u8) -> i64 {
    -((10f64.powf(-(phred as f64) / 10.0) * ERROR_UNIT).round() as i64)
}

static PHRED_SCALE: LazyLock<Scale> = LazyLock::new(|| Scale {
    score: std::array::from_fn(|c| c as i64),
    averaging: Averaging::Phred,
});

static ERROR_SCALE: LazyLock<Scale> = LazyLock::new(|| Scale {
    score: std::array::from_fn(|c| error_score((c as u8).saturating_sub(PHRED_OFFSET as u8))),
    averaging: Averaging::ErrorProbability,
});

impl Scale {
    pub fn get(averaging: Averaging) -> &'static Scale {
        match averaging {
            Averaging::Phred => &PHRED_SCALE,
            Averaging::ErrorProbability => &ERROR_SCALE,
        }
    }

    /// Score of a quality character.
    #[inline]
    pub fn score(&self, q: u8) -> i64 {
        self.score[q as usize]
    }

    /// Score of a base with phred quality `phred`.
    pub fn threshold(&self, phred: u8) -> i64 {
        match self.averaging {
            Averaging::Phred => PHRED_OFFSET + phred as i64,
            Averaging::ErrorProbability => error_score(phred),
        }
    }
}

/// Mean quality of a read as phred, averaging error probabilities. `None` for an empty read.
pub fn mean_error_quality(qual: &[u8]) -> Option<f64> {
    if qual.is_empty() {
        return None;
    }
    Some(phred_of_error(error_sum(qual) / qual.len() as f64))
}

/// Sum of the error probabilities of `qual`.
pub fn error_sum(qual: &[u8]) -> f64 {
    -qual.iter().map(|&q| ERROR_SCALE.score(q)).sum::<i64>() as f64 / ERROR_UNIT
}

/// Phred value of an error probability (capped at 60 for error-free input).
pub fn phred_of_error(p: f64) -> f64 {
    if p <= 1e-6 { 60.0 } else { -10.0 * p.log10() }
}

const PHRED_OFFSET: i64 = 33;

/// Applies fixed trimming (`front`/`tail` bases) and then the enabled sliding-window cuts,
/// with window means taken on `scale`.
///
/// Returns the kept region `[start, end)` of the read, or `None` if nothing is left.
///
/// Differences to fastp (deliberate bug fixes):
/// - when fixed trimming is combined with `cut_front`/`cut_tail`, fastp drops `window - 1`
///   extra good bases even if the first window already passes; here nothing extra is cut;
/// - `cut_right` never keeps bases from the fixed-trimmed tail.
pub fn trim_and_cut(
    seq: &[u8],
    qual: &[u8],
    front: usize,
    tail: usize,
    cut: &QualityCut,
    scale: &Scale,
) -> Option<(usize, usize)> {
    if front == 0 && tail == 0 && !cut.is_enabled() {
        return Some((0, seq.len()));
    }
    let l = seq.len() as isize;
    let tail = tail as isize;
    let mut front = front as isize;
    let mut rlen = l - front - tail;
    if rlen < 0 {
        return None;
    }
    if !cut.is_enabled() {
        return Some((front as usize, (front + rlen) as usize));
    }

    let q = |i: isize| scale.score(qual[i as usize]);

    if let Some(w) = cut.front {
        let ws = w.size as isize;
        if l - front - tail - ws <= 0 {
            return None;
        }
        let threshold = ws as i64 * scale.threshold(w.quality);
        let mut total: i64 = (0..ws - 1).map(|i| q(front + i)).sum();
        let mut s = front;
        while s + ws < l - tail {
            total += q(s + ws - 1);
            if s > front {
                total -= q(s - 1);
            }
            if total >= threshold {
                break;
            }
            s += 1;
        }
        // drop every base covered by a failing window
        if s > front {
            s = s + ws - 1;
        }
        while s < l && seq[s as usize] == b'N' {
            s += 1;
        }
        front = s;
        rlen = l - front - tail;
    }

    if let Some(w) = cut.right {
        let ws = w.size as isize;
        if l - front - tail - ws <= 0 {
            return None;
        }
        let base_threshold = scale.threshold(w.quality);
        let threshold = ws as i64 * base_threshold;
        let mut total: i64 = (0..ws - 1).map(|i| q(front + i)).sum();
        let mut s = front;
        let mut found = false;
        while s + ws < l - tail {
            total += q(s + ws - 1);
            if s > front {
                total -= q(s - 1);
            }
            if total < threshold {
                found = true;
                break;
            }
            s += 1;
        }
        if found {
            // keep the good bases at the start of the failing window
            while s < l - 1 && q(s) >= base_threshold {
                s += 1;
            }
            rlen = (s - front).min(l - front - tail);
        }
    } else if let Some(w) = cut.tail {
        let ws = w.size as isize;
        if l - front - tail - ws <= 0 {
            return None;
        }
        let threshold = ws as i64 * scale.threshold(w.quality);
        let mut t = l - tail - 1;
        let mut total: i64 = (0..ws - 1).map(|i| q(t - i)).sum();
        while t - ws >= front {
            total += q(t - ws + 1);
            if t < l - tail - 1 {
                total -= q(t + 1);
            }
            if total >= threshold {
                break;
            }
            t -= 1;
        }
        if t < l - tail - 1 {
            t = t - ws + 1;
        }
        while t >= 0 && seq[t as usize] == b'N' {
            t -= 1;
        }
        rlen = t - front + 1;
    }

    if rlen <= 0 || front >= l - 1 {
        return None;
    }
    Some((front as usize, (front + rlen) as usize))
}

/// Low-quality stretches of a read: every window of `window` bases whose mean quality is
/// below `quality` is marked, overlapping or touching windows are merged, and each stretch is
/// then shrunk to start and end at a base below `quality`, so that good bases at the edge of a
/// window are kept. A read shorter than the window is judged as one window.
///
/// Returns the stretches as `[start, end)` ranges in read order.
pub fn low_quality_regions(
    qual: &[u8],
    window: usize,
    quality: u8,
    scale: &Scale,
) -> Vec<(usize, usize)> {
    let l = qual.len();
    if l == 0 || window == 0 {
        return Vec::new();
    }
    let w = window.min(l);
    let base_threshold = scale.threshold(quality);
    let threshold = w as i64 * base_threshold;
    let score = |i: usize| scale.score(qual[i]);
    let mut total: i64 = (0..w).map(score).sum();
    let mut marked: Vec<(usize, usize)> = Vec::new();
    for s in 0..=l - w {
        if s > 0 {
            total += score(s + w - 1) - score(s - 1);
        }
        if total < threshold {
            match marked.last_mut() {
                Some(last) if s <= last.1 => last.1 = s + w,
                _ => marked.push((s, s + w)),
            }
        }
    }
    marked
        .into_iter()
        .filter_map(|(mut s, mut e)| {
            while s < e && score(s) >= base_threshold {
                s += 1;
            }
            while e > s && score(e - 1) >= base_threshold {
                e -= 1;
            }
            (s < e).then_some((s, e))
        })
        .collect()
}

/// The parts of `[0, len)` outside the given sorted, non-overlapping `regions`.
pub fn complement(len: usize, regions: &[(usize, usize)]) -> Vec<(usize, usize)> {
    let mut out = Vec::with_capacity(regions.len() + 1);
    let mut pos = 0;
    for &(s, e) in regions {
        if s > pos {
            out.push((pos, s));
        }
        pos = pos.max(e);
    }
    if pos < len {
        out.push((pos, len));
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn phred() -> &'static Scale {
        Scale::get(Averaging::Phred)
    }

    fn errp() -> &'static Scale {
        Scale::get(Averaging::ErrorProbability)
    }

    const W4Q20: Window = Window {
        size: 4,
        quality: 20,
    };

    #[test]
    fn no_trimming_returns_whole_read() {
        let cut = QualityCut::default();
        assert_eq!(
            trim_and_cut(b"ACGT", b"IIII", 0, 0, &cut, phred()),
            Some((0, 4))
        );
    }

    #[test]
    fn fixed_trimming() {
        let cut = QualityCut::default();
        assert_eq!(
            trim_and_cut(b"ACGTACGT", b"IIIIIIII", 2, 1, &cut, phred()),
            Some((2, 7))
        );
        assert_eq!(trim_and_cut(b"ACGT", b"IIII", 3, 2, &cut, phred()), None);
        // trimming everything leaves an empty read, which the length filter removes
        assert_eq!(
            trim_and_cut(b"ACGT", b"IIII", 2, 2, &cut, phred()),
            Some((2, 2))
        );
    }

    /// Test vector from fastp's `Filter::test`.
    #[test]
    fn fastp_front_and_tail_cut() {
        let seq = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTT";
        let qual = b"/////CCCCCCCCCCCC////CCCCCCCCCCCCCC////E";
        let cut = QualityCut {
            front: Some(W4Q20),
            tail: Some(W4Q20),
            right: None,
        };
        let (s, e) = trim_and_cut(seq, qual, 0, 1, &cut, phred()).unwrap();
        assert_eq!(&seq[s..e], b"CCCCCCCCCCCCCCCCCCCCCCCCCCCC");
        assert_eq!(&qual[s..e], b"CCCCCCCCCCC////CCCCCCCCCCCCC");
    }

    #[test]
    fn cut_right_drops_everything_after_bad_window() {
        let seq = b"ACGTACGTACGTACGTACGT";
        let qual = b"IIIIIIIIII#####IIIII";
        let cut = QualityCut {
            front: None,
            tail: None,
            right: Some(W4Q20),
        };
        let (s, e) = trim_and_cut(seq, qual, 0, 0, &cut, phred()).unwrap();
        assert_eq!((s, e), (0, 10));
    }

    #[test]
    fn good_first_window_with_fixed_trim_cuts_nothing_extra() {
        let seq = b"ACGTACGTACGTACGTACGT";
        let qual = b"IIIIIIIIIIIIIIIIIIII";
        let cut = QualityCut {
            front: Some(W4Q20),
            tail: Some(W4Q20),
            right: None,
        };
        assert_eq!(trim_and_cut(seq, qual, 2, 2, &cut, phred()), Some((2, 18)));
    }

    #[test]
    fn all_low_quality_is_removed() {
        let seq = b"ACGTACGTACGT";
        let qual = b"############";
        let cut = QualityCut {
            front: Some(W4Q20),
            tail: None,
            right: None,
        };
        assert_eq!(trim_and_cut(seq, qual, 0, 0, &cut, phred()), None);
    }

    #[test]
    fn error_probability_means() {
        assert!((mean_error_quality(b"++++").unwrap() - 10.0).abs() < 1e-6);
        // Q30 and Q10: plain mean 20, error-probability mean about Q13
        let q = mean_error_quality(b"?+").unwrap();
        assert!((q - 12.967).abs() < 0.01, "{q}");
        assert_eq!(mean_error_quality(b""), None);
        // a window of exactly the threshold quality passes
        let s = errp();
        assert!(4 * s.score(b'5') >= 4 * s.threshold(20));
        assert!(s.score(b'4') < s.threshold(20));
    }

    #[test]
    fn error_probability_cut_is_stricter() {
        // one Q2 base among Q40: the plain mean of a 4-base window stays above Q20, the
        // error-probability mean drops to about Q11
        let seq = b"ACGTACGTACGT";
        let cut = QualityCut {
            front: Some(W4Q20),
            tail: None,
            right: None,
        };
        let qual = b"II#IIIIIIIII";
        assert_eq!(trim_and_cut(seq, qual, 0, 0, &cut, phred()), Some((0, 12)));
        assert_eq!(trim_and_cut(seq, qual, 0, 0, &cut, errp()), Some((6, 12)));
    }

    #[test]
    fn low_quality_stretches() {
        // two dips in a good read; fastplong's detection loses the second one (it marks
        // everything after the first dip)
        let mut qual = vec![b'?'; 1000];
        for q in &mut qual[200..260] {
            *q = b'#';
        }
        for q in &mut qual[600..660] {
            *q = b'#';
        }
        let r = low_quality_regions(&qual, 50, 10, phred());
        assert_eq!(r, vec![(200, 260), (600, 660)]);
        assert_eq!(
            complement(1000, &r),
            vec![(0, 200), (260, 600), (660, 1000)]
        );
        assert_eq!(low_quality_regions(&qual, 50, 10, errp()), r);
    }

    #[test]
    fn low_quality_edges_and_short_reads() {
        assert!(low_quality_regions(&[b'?'; 300], 100, 10, errp()).is_empty());
        // a few bad bases spread over a window pull the error-probability mean down
        let mut qual = vec![b'5'; 300]; // Q20
        for i in (100..140).step_by(4) {
            qual[i] = b'!'; // Q0
        }
        assert!(low_quality_regions(&qual, 100, 10, phred()).is_empty());
        assert_eq!(
            low_quality_regions(&qual, 100, 10, errp()),
            vec![(100, 137)]
        );
        // a read shorter than the window is judged as a whole
        assert_eq!(
            low_quality_regions(b"##########", 100, 10, errp()),
            vec![(0, 10)]
        );
        assert!(low_quality_regions(b"??????????", 100, 10, errp()).is_empty());
        assert_eq!(complement(10, &[(0, 10)]), Vec::<(usize, usize)>::new());
    }
}
