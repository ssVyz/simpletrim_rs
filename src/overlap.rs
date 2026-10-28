//! Read-pair overlap analysis (port of fastp's `OverlapAnalysis`) and overlap-based
//! adapter trimming.

use crate::matcher::{
    count_mismatches, count_mismatches_bounded, diff_with_one_insertion, reverse_complement_into,
};

/// Thresholds for accepting an overlap between read 1 and the reverse complement of read 2.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OverlapParams {
    /// Minimum overlap length (fastp `--overlap_len_require`, 30).
    pub min_len: usize,
    /// Maximum mismatches in the overlap (fastp `--overlap_diff_limit`, 5).
    pub max_mismatches: usize,
    /// Maximum mismatch fraction of the overlap (fastp `--overlap_diff_percent_limit`, 20%).
    pub max_error_rate: f64,
    /// Also try alignments with one insertion/deletion (fastp `--allow_gap_overlap_trimming`).
    pub allow_gap: bool,
}

impl Default for OverlapParams {
    fn default() -> Self {
        OverlapParams {
            min_len: 30,
            max_mismatches: 5,
            max_error_rate: 0.2,
            allow_gap: false,
        }
    }
}

/// Result of an overlap search.
///
/// `offset` is the position of reverse-complemented read 2 relative to read 1; a negative
/// offset means the insert is shorter than the reads, i.e. both reads run into adapter.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub struct Overlap {
    pub overlapped: bool,
    pub offset: isize,
    pub overlap_len: usize,
    pub diff: usize,
    pub has_gap: bool,
}

/// Only the first bases of an overlap must stay within the mismatch limit; the tail of long
/// overlaps is usually low quality (fastp `complete_compare_require`).
const COMPLETE_COMPARE_REQUIRE: usize = 50;

fn accept_no_gap(a: &[u8], b: &[u8], len: usize, limit: usize) -> Option<usize> {
    let protected = len.min(COMPLETE_COMPARE_REQUIRE);
    let diff = count_mismatches_bounded(&a[..protected], &b[..protected], limit);
    if diff > limit {
        return None;
    }
    if len > COMPLETE_COMPARE_REQUIRE {
        return Some(count_mismatches(&a[..len], &b[..len]));
    }
    Some(diff)
}

fn gapped_diff(a: &[u8], b: &[u8], overlap_len: usize, limit: usize) -> Option<usize> {
    let cmplen = overlap_len.saturating_sub(1);
    diff_with_one_insertion(a, b, cmplen, limit as isize)
        .or_else(|| diff_with_one_insertion(b, a, cmplen, limit as isize))
}

/// Finds the overlap between `r1` and the reverse complement of `r2`.
/// `rc_buf` is scratch space for the reverse complement.
pub fn analyze(r1: &[u8], r2: &[u8], p: &OverlapParams, rc_buf: &mut Vec<u8>) -> Overlap {
    reverse_complement_into(r2, rc_buf);
    let s1 = r1;
    let s2 = rc_buf.as_slice();
    let len1 = s1.len() as isize;
    let len2 = s2.len() as isize;
    let req = p.min_len as isize;
    let limit_for = |overlap_len: usize| {
        p.max_mismatches
            .min((overlap_len as f64 * p.max_error_rate) as usize)
    };
    let found = |offset: isize, overlap_len: usize, diff: usize, has_gap: bool| Overlap {
        overlapped: true,
        offset,
        overlap_len,
        diff,
        has_gap,
    };

    // read 2 shifted right of read 1: insert longer than the reads, no adapter
    let mut offset = 0;
    while offset < len1 - req {
        let overlap_len = (len1 - offset).min(len2) as usize;
        let o = offset as usize;
        if let Some(diff) = accept_no_gap(&s1[o..], s2, overlap_len, limit_for(overlap_len)) {
            return found(offset, overlap_len, diff, false);
        }
        offset += 1;
    }

    // read 2 shifted left: insert shorter than the reads, adapter sequenced
    let mut offset = 0;
    while offset > -(len2 - req) {
        let overlap_len = len1.min(len2 - offset.abs()) as usize;
        let o = (-offset) as usize;
        if let Some(diff) = accept_no_gap(s1, &s2[o..], overlap_len, limit_for(overlap_len)) {
            return found(offset, overlap_len, diff, false);
        }
        offset -= 1;
    }

    if p.allow_gap {
        let mut offset = 0;
        while offset < len1 - req {
            let overlap_len = (len1 - offset).min(len2) as usize;
            let o = offset as usize;
            if let Some(diff) = gapped_diff(&s1[o..], s2, overlap_len, limit_for(overlap_len)) {
                return found(offset, overlap_len, diff, true);
            }
            offset += 1;
        }
        let mut offset = 0;
        while offset > -(len2 - req) {
            let overlap_len = len1.min(len2 - offset.abs()) as usize;
            let o = (-offset) as usize;
            if let Some(diff) = gapped_diff(s1, &s2[o..], overlap_len, limit_for(overlap_len)) {
                return found(offset, overlap_len, diff, true);
            }
            offset -= 1;
        }
    }

    Overlap::default()
}

/// Lengths to keep for both reads when the overlap shows adapter read-through
/// (fastp `trimByOverlapAnalysis`). `front1`/`front2` are the bases already removed from the
/// 5' ends, which are added back to the mate's insert length.
pub fn adapter_cut(
    len1: usize,
    len2: usize,
    ov: &Overlap,
    front1: usize,
    front2: usize,
) -> Option<(usize, usize)> {
    if ov.overlapped && ov.offset < 0 {
        Some((
            len1.min(ov.overlap_len + front2),
            len2.min(ov.overlap_len + front1),
        ))
    } else {
        None
    }
}

/// Insert size implied by an overlap (fastp `statInsertSize`), or `None` if not overlapped.
pub fn insert_size(len1: usize, len2: usize, ov: &Overlap) -> Option<usize> {
    if !ov.overlapped {
        return None;
    }
    Some(if ov.offset > 0 {
        len1 + len2 - ov.overlap_len
    } else {
        ov.overlap_len
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    fn rc(s: &[u8]) -> Vec<u8> {
        let mut out = Vec::new();
        reverse_complement_into(s, &mut out);
        out
    }

    /// Test vectors from fastp's `OverlapAnalysis::test`.
    #[test]
    fn fastp_overlap() {
        let r1 = b"CAGCGCCTACGGGCCCCTTTTTCTGCGCGACCGCGTGGCTGTGGGCGCGGATGCCTTTGAGCGCGGTGACTTCTCACTGCGTATCGAGC";
        let r2 = b"ACCTCCAGCGGCTCGATACGCAGTGAGAAGTCACCGCGCTCAAAGGCATCCGCGCCCACAGCCACGCGGTCGCGCAGAAAAAGGGGTCC";
        let p = OverlapParams {
            min_len: 30,
            max_mismatches: 2,
            max_error_rate: 0.2,
            allow_gap: false,
        };
        let ov = analyze(r1, r2, &p, &mut Vec::new());
        assert!(ov.overlapped);
        assert_eq!((ov.offset, ov.overlap_len, ov.diff), (10, 79, 1));
    }

    #[test]
    fn fastp_late_mismatches_are_accepted() {
        let r1 = [vec![b'A'; 50], vec![b'C'; 30]].concat();
        let rcr2 = [vec![b'A'; 50], vec![b'G'; 30]].concat();
        let r2 = rc(&rcr2);
        let p = OverlapParams {
            min_len: 30,
            max_mismatches: 0,
            max_error_rate: 0.0,
            allow_gap: false,
        };
        let ov = analyze(&r1, &r2, &p, &mut Vec::new());
        assert!(ov.overlapped);
        assert_eq!((ov.offset, ov.overlap_len, ov.diff), (0, 80, 30));
    }

    #[test]
    fn short_insert_trims_adapters() {
        let insert = b"GATTACAGATTACACCATGGTTGACCATGGATCCAAGTTCGGATCAGT"; // 48 bp
        let a1 = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
        let a2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTA";
        let r1 = [insert.as_slice(), &a1[..30]].concat();
        let r2 = [rc(insert).as_slice(), &a2[..30]].concat();
        let ov = analyze(&r1, &r2, &OverlapParams::default(), &mut Vec::new());
        assert!(ov.overlapped && ov.offset < 0);
        assert_eq!(adapter_cut(r1.len(), r2.len(), &ov, 0, 0), Some((48, 48)));
        assert_eq!(insert_size(r1.len(), r2.len(), &ov), Some(48));
    }

    #[test]
    fn unrelated_reads_do_not_overlap() {
        let r1 = b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGTTGCAACGT";
        let r2 = b"GGGCCCAAATTTGGGCCCAAATTTGGGCCCAAATTTGGGCCCAA";
        let ov = analyze(r1, r2, &OverlapParams::default(), &mut Vec::new());
        assert!(!ov.overlapped);
        assert_eq!(adapter_cut(r1.len(), r2.len(), &ov, 0, 0), None);
    }

    #[test]
    fn gapped_overlap() {
        let insert = b"GATTACAGATTACACCATGGTTGACCATGGATCCAAGTTCGGATCAGTTGCA"; // 52 bp
        let a1 = b"AGATCGGAAGAGCACACGTCTGAACTCC";
        let a2 = b"AGATCGGAAGAGCGTCGTGTAGGGAAAG";
        let r1 = [insert.as_slice(), a1.as_slice()].concat();
        // read 2 with one base deleted in the middle of the insert
        let mut ins2 = rc(insert);
        ins2.remove(20);
        let r2 = [ins2.as_slice(), a2.as_slice()].concat();
        let strict = OverlapParams {
            max_mismatches: 2,
            ..OverlapParams::default()
        };
        let no_gap = analyze(&r1, &r2, &strict, &mut Vec::new());
        assert!(!no_gap.overlapped || no_gap.offset >= 0);
        let gap = OverlapParams {
            allow_gap: true,
            ..strict
        };
        let ov = analyze(&r1, &r2, &gap, &mut Vec::new());
        assert!(ov.overlapped && ov.has_gap && ov.offset < 0);
    }
}
