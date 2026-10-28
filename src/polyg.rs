//! Poly-G tail trimming (port of fastp's `PolyX::trimPolyG`).

const ALLOW_ONE_MISMATCH_FOR_EACH: usize = 8;
const MAX_MISMATCH: usize = 5;

/// Detects a poly-G tail of at least `min_len` bases (one mismatch allowed per 8 bases,
/// at most 5). Returns the read length after removing the tail, or `None` if there is none.
pub fn poly_g_trim_len(seq: &[u8], min_len: usize) -> Option<usize> {
    let rlen = seq.len();
    if rlen == 0 || min_len == 0 {
        return None;
    }
    let mut mismatch = 0;
    let mut first_g = rlen - 1;
    let mut i = 0;
    while i < rlen {
        let pos = rlen - i - 1;
        if seq[pos] != b'G' {
            mismatch += 1;
        } else {
            first_g = pos;
        }
        let allowed = (i + 1) / ALLOW_ONE_MISMATCH_FOR_EACH;
        if mismatch > MAX_MISMATCH || (mismatch > allowed && i + 1 >= min_len) {
            break;
        }
        i += 1;
    }
    (i >= min_len).then_some(first_g)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn trims_poly_g_tail() {
        let seq = b"ACATACATACATGGGGGGGGGGGGGG";
        assert_eq!(poly_g_trim_len(seq, 10), Some(12));
    }

    #[test]
    fn tolerates_a_mismatch() {
        let seq = b"ACATACATACATGGGGGGAGGGGGGGGG";
        assert_eq!(poly_g_trim_len(seq, 10), Some(12));
    }

    #[test]
    fn extends_through_mismatch_onto_g_like_fastp() {
        // the G at position 10 is within the mismatch budget and becomes the cut point
        let seq = b"ACGTACGTACGTGGGGGGGGGGGGGG";
        assert_eq!(poly_g_trim_len(seq, 10), Some(10));
    }

    #[test]
    fn short_tail_is_kept() {
        assert_eq!(poly_g_trim_len(b"ACGTACGTACGTACGGGGG", 10), None);
        assert_eq!(poly_g_trim_len(b"", 10), None);
    }

    #[test]
    fn whole_read_poly_g() {
        assert_eq!(poly_g_trim_len(b"GGGGGGGGGGGG", 10), Some(0));
    }
}
