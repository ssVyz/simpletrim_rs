//! Poly-X tail trimming (port of fastp's `PolyX::trimPolyX`, also used by fastplong).

const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 8;
const MAX_MISMATCH: usize = 5;

/// Index of a base in the A, T, C, G counters; `N` counts for every base.
fn base_index(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'T' | b'U' => Some(1),
        b'C' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Detects a homopolymer tail (any base) of at least `min_len` bases at the 3' end, allowing
/// one mismatch per 8 bases (at most 5). Returns the read length after removing the tail, or
/// `None` if there is none.
///
/// Difference to fastp: a read that is a homopolymer from end to end is trimmed completely;
/// fastp reads before the start of the read in that case.
pub fn poly_x_trim_len(seq: &[u8], min_len: usize) -> Option<usize> {
    let rlen = seq.len();
    if rlen == 0 || min_len == 0 {
        return None;
    }
    let mut counts = [0usize; 4];
    let mut pos = 0;
    while pos < rlen {
        match seq[rlen - pos - 1] {
            b'N' => counts.iter_mut().for_each(|c| *c += 1),
            b => {
                if let Some(i) = base_index(b) {
                    counts[i] += 1;
                }
            }
        }
        let cmp = pos + 1;
        let allowed = MAX_MISMATCH.min(cmp / ALLOW_ONE_MISMATCH_FOR_EACH);
        let need_to_break = counts.iter().all(|&c| cmp - c > allowed);
        if need_to_break && (pos >= ALLOW_ONE_MISMATCH_FOR_EACH || pos + 2 >= min_len) {
            break;
        }
        pos += 1;
    }
    if pos + 1 < min_len {
        return None;
    }
    // the most frequent base (first in A, T, C, G order on ties)
    let poly = (0..4).fold(0, |best, b| if counts[b] > counts[best] { b } else { best });
    // the tail ends at the last occurrence of the poly base within the scanned part
    let mut p = pos.min(rlen - 1) as isize;
    while p >= 0 && base_index(seq[rlen - p as usize - 1]) != Some(poly) {
        p -= 1;
    }
    let new_len = (rlen as isize - p - 1) as usize;
    (new_len < rlen).then_some(new_len)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test vector from fastp's `PolyX::test`.
    #[test]
    fn fastp_vector() {
        let seq = b"ATTTTAAAAAAAAAATAAAAAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAT";
        assert_eq!(poly_x_trim_len(seq, 10), Some(5));
    }

    #[test]
    fn tails_of_each_base() {
        assert_eq!(poly_x_trim_len(b"ACGTACGTACGTCCCCCCCCCCCC", 10), Some(12));
        assert_eq!(poly_x_trim_len(b"ACGTACGTACGTTTTTTTTTTTTT", 10), Some(11));
        // too short
        assert_eq!(poly_x_trim_len(b"ACGTACGTACGTAAAAAAA", 10), None);
        assert_eq!(poly_x_trim_len(b"ACGTACGTACGTACGTACGT", 10), None);
    }

    #[test]
    fn whole_read_homopolymer() {
        assert_eq!(poly_x_trim_len(b"AAAAAAAAAAAAAAAA", 10), Some(0));
        assert_eq!(poly_x_trim_len(b"", 10), None);
    }
}
