//! Fixed front/tail trimming and sliding-window quality cutting (port of fastp's `trimAndCut`).

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

const PHRED_OFFSET: i64 = 33;

/// Applies fixed trimming (`front`/`tail` bases) and then the enabled sliding-window cuts.
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

    let q = |i: isize| qual[i as usize] as i64;

    if let Some(w) = cut.front {
        let ws = w.size as isize;
        if l - front - tail - ws <= 0 {
            return None;
        }
        let threshold = ws as i64 * (PHRED_OFFSET + w.quality as i64);
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
        let base_threshold = PHRED_OFFSET + w.quality as i64;
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
        let threshold = ws as i64 * (PHRED_OFFSET + w.quality as i64);
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

#[cfg(test)]
mod tests {
    use super::*;

    const W4Q20: Window = Window {
        size: 4,
        quality: 20,
    };

    #[test]
    fn no_trimming_returns_whole_read() {
        let cut = QualityCut::default();
        assert_eq!(trim_and_cut(b"ACGT", b"IIII", 0, 0, &cut), Some((0, 4)));
    }

    #[test]
    fn fixed_trimming() {
        let cut = QualityCut::default();
        assert_eq!(
            trim_and_cut(b"ACGTACGT", b"IIIIIIII", 2, 1, &cut),
            Some((2, 7))
        );
        assert_eq!(trim_and_cut(b"ACGT", b"IIII", 3, 2, &cut), None);
        // trimming everything leaves an empty read, which the length filter removes
        assert_eq!(trim_and_cut(b"ACGT", b"IIII", 2, 2, &cut), Some((2, 2)));
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
        let (s, e) = trim_and_cut(seq, qual, 0, 1, &cut).unwrap();
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
        let (s, e) = trim_and_cut(seq, qual, 0, 0, &cut).unwrap();
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
        assert_eq!(trim_and_cut(seq, qual, 2, 2, &cut), Some((2, 18)));
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
        assert_eq!(trim_and_cut(seq, qual, 0, 0, &cut), None);
    }
}
