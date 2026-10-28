//! Sequence comparison primitives shared by adapter and overlap detection.

use std::cell::RefCell;

/// Complement of a nucleotide; anything outside `ACGT` (any case) becomes `N`.
pub fn complement(base: u8) -> u8 {
    match base {
        b'A' | b'a' => b'T',
        b'T' | b't' => b'A',
        b'C' | b'c' => b'G',
        b'G' | b'g' => b'C',
        _ => b'N',
    }
}

/// Writes the reverse complement of `src` into `dst` (replacing its contents).
pub fn reverse_complement_into(src: &[u8], dst: &mut Vec<u8>) {
    dst.clear();
    dst.extend(src.iter().rev().map(|&b| complement(b)));
}

/// Number of differing positions over the common length of `a` and `b`.
pub fn count_mismatches(a: &[u8], b: &[u8]) -> usize {
    a.iter().zip(b).filter(|(x, y)| x != y).count()
}

/// Like [`count_mismatches`], but stops once `limit` is exceeded (the result is then `limit + 1`).
pub fn count_mismatches_bounded(a: &[u8], b: &[u8], limit: usize) -> usize {
    let mut n = 0;
    for (x, y) in a.iter().zip(b) {
        if x != y {
            n += 1;
            if n > limit {
                break;
            }
        }
    }
    n
}

thread_local! {
    static SCRATCH: RefCell<(Vec<usize>, Vec<usize>)> = const { RefCell::new((Vec::new(), Vec::new())) };
}

/// Compares `ins` (`cmplen + 1` bases, containing one extra base) against `normal`
/// (`cmplen` bases), allowing the extra base at any position `1..cmplen`.
///
/// Port of fastp's `Matcher::matchWithOneInsertion` / `diffWithOneInsertion`, including their
/// early-termination rules. With `first_match` the first placement within `limit` is returned,
/// otherwise the smallest mismatch count. Returns `None` if no placement is within `limit`.
fn one_insertion(
    ins: &[u8],
    normal: &[u8],
    cmplen: usize,
    limit: isize,
    first_match: bool,
) -> Option<usize> {
    if limit < 0 || cmplen < 2 || ins.len() <= cmplen || normal.len() < cmplen {
        return None;
    }
    let limit = limit as usize;
    SCRATCH.with(|cell| {
        let (left, right) = &mut *cell.borrow_mut();
        left.resize(cmplen, 0);
        right.resize(cmplen, 0);

        // left[i]: mismatches of ins[0..=i] vs normal[0..=i]
        // right[i]: mismatches of ins[i+1..=cmplen] vs normal[i..cmplen]
        let right_last = (ins[cmplen] != normal[cmplen - 1]) as usize;
        left[0] = (ins[0] != normal[0]) as usize;
        for i in 1..cmplen {
            left[i] = left[i - 1] + (ins[i] != normal[i]) as usize;
            if left[i] + right_last > limit {
                // later entries are never read: the final loop stops at index i
                break;
            }
        }
        right[cmplen - 1] = right_last;
        let mut right_valid_from = cmplen - 1;
        for i in (0..cmplen - 1).rev() {
            right[i] = right[i + 1] + (ins[i + 1] != normal[i]) as usize;
            right_valid_from = i;
            if right[i] + left[0] > limit {
                break;
            }
        }

        let mut best: Option<usize> = None;
        for i in 1..cmplen {
            if left[i - 1] + right_last > limit {
                break;
            }
            let r = if i < right_valid_from {
                limit + 1
            } else {
                right[i]
            };
            let diff = left[i - 1] + r;
            if diff <= limit {
                if first_match {
                    return Some(diff);
                }
                best = Some(best.map_or(diff, |b| b.min(diff)));
            }
        }
        best
    })
}

/// True if `ins` matches `normal` over `cmplen` bases with one inserted base in `ins`
/// and at most `limit` mismatches.
pub fn matches_with_one_insertion(ins: &[u8], normal: &[u8], cmplen: usize, limit: isize) -> bool {
    one_insertion(ins, normal, cmplen, limit, true).is_some()
}

/// Smallest mismatch count of `ins` vs `normal` with one inserted base in `ins`, if within `limit`.
pub fn diff_with_one_insertion(
    ins: &[u8],
    normal: &[u8],
    cmplen: usize,
    limit: isize,
) -> Option<usize> {
    one_insertion(ins, normal, cmplen, limit, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn reverse_complement_works() {
        let mut out = Vec::new();
        reverse_complement_into(b"AACGTN", &mut out);
        assert_eq!(out, b"NACGTT");
    }

    #[test]
    fn bounded_mismatch_stops_early() {
        assert_eq!(count_mismatches(b"AAAA", b"ATTA"), 2);
        assert_eq!(count_mismatches_bounded(b"AAAA", b"TTTT", 1), 2);
        assert_eq!(count_mismatches_bounded(b"AAAA", b"AAAT", 1), 1);
    }

    #[test]
    fn one_insertion_matches() {
        // "ACGTXACGT" vs "ACGTACGT": the X is the inserted base
        let ins = b"ACGTTACGTA";
        let normal = b"ACGTACGTA";
        assert!(matches_with_one_insertion(ins, normal, 9, 0));
        assert_eq!(diff_with_one_insertion(ins, normal, 9, 1), Some(0));
        // unrelated sequences do not match
        assert!(!matches_with_one_insertion(
            b"GGGGGGGGGG",
            b"ACACACACA",
            9,
            1
        ));
        // negative limits never match
        assert!(!matches_with_one_insertion(ins, normal, 9, -1));
    }

    #[test]
    fn one_insertion_reports_minimum() {
        // insertion near the start, one mismatch late in the sequence
        let ins = b"AGCCCCCCCCCCT";
        let normal = b"ACCCCCCCCCCA";
        assert_eq!(diff_with_one_insertion(ins, normal, 12, 2), Some(1));
        assert_eq!(diff_with_one_insertion(ins, normal, 12, 0), None);
    }
}
