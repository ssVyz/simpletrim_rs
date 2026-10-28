//! Sequence-based adapter trimming (port of fastp's `AdapterTrimmer::trimBySequence`).

use crate::matcher::{count_mismatches_bounded, matches_with_one_insertion};

/// Tolerances for matching a known adapter sequence against a read.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MatchParams {
    /// Minimum number of adapter bases that must overlap the read end (fastp: 4).
    /// `None` uses fastp's rule: 4, raised to 5/6 for adapter lists longer than 16/256 entries.
    pub min_overlap: Option<usize>,
    /// Allowed mismatches per compared base (fastp: one per 8 bases = 0.125).
    pub error_rate: f64,
    /// Absolute cap on mismatches, independent of the compared length.
    pub max_mismatches: Option<usize>,
}

impl Default for MatchParams {
    fn default() -> Self {
        MatchParams {
            min_overlap: None,
            error_rate: 0.125,
            max_mismatches: None,
        }
    }
}

impl MatchParams {
    /// Allowed mismatches when comparing `cmplen` bases.
    fn allowed(&self, cmplen: usize) -> isize {
        let by_rate = (cmplen as f64 * self.error_rate + 1e-9).floor() as isize;
        match self.max_mismatches {
            Some(cap) => by_rate.min(cap as isize),
            None => by_rate,
        }
    }

    fn min_overlap_for(&self, list_len: usize) -> usize {
        self.min_overlap.unwrap_or(match list_len {
            n if n > 256 => 6,
            n if n > 16 => 5,
            _ => 4,
        })
    }
}

/// Searches `adapter` in `seq` and returns the length to keep, or `None` if not found.
///
/// The search allows the adapter to start slightly before the read (Illumina adapter dimers
/// often lack the first A due to A-tailing), then tries a single insertion and a single
/// deletion. Unlike fastp, the gapped searches compare at the candidate position rather than
/// at the start of the read.
pub fn find_adapter(
    seq: &[u8],
    adapter: &[u8],
    min_overlap: usize,
    p: &MatchParams,
) -> Option<usize> {
    let rlen = seq.len() as isize;
    let alen = adapter.len() as isize;
    let mreq = min_overlap as isize;
    if alen < mreq || alen == 0 {
        return None;
    }

    let start: isize = match alen {
        16.. => -4,
        12.. => -3,
        8.. => -2,
        _ => 0,
    };
    let mut pos = start;
    while pos < rlen - mreq {
        let cmplen = (rlen - pos).min(alen);
        let allowed = p.allowed(cmplen as usize);
        let skip = (-pos).max(0);
        if allowed >= 0 && skip < cmplen {
            let a = &adapter[skip as usize..cmplen as usize];
            let r = &seq[(skip + pos) as usize..(pos + cmplen) as usize];
            if count_mismatches_bounded(a, r, allowed as usize) <= allowed as usize {
                return Some(pos.max(0) as usize);
            }
        }
        pos += 1;
    }

    // one inserted base in the read
    let mut pos = 0;
    while pos < rlen - mreq - 1 {
        let cmplen = (rlen - pos - 1).min(alen);
        let allowed = p.allowed(cmplen as usize) - 1;
        if matches_with_one_insertion(&seq[pos as usize..], adapter, cmplen as usize, allowed) {
            return Some(pos as usize);
        }
        pos += 1;
    }

    // one deleted base in the read
    let mut pos = 0;
    while pos < rlen - mreq {
        let cmplen = (rlen - pos).min(alen - 1);
        let allowed = p.allowed(cmplen as usize) - 1;
        if matches_with_one_insertion(adapter, &seq[pos as usize..], cmplen as usize, allowed) {
            return Some(pos as usize);
        }
        pos += 1;
    }
    None
}

/// Trims a single configured adapter; returns the new length if the adapter was found.
pub fn trim_single(seq: &[u8], adapter: &[u8], p: &MatchParams) -> Option<usize> {
    find_adapter(seq, adapter, p.min_overlap_for(1), p)
}

/// Trims every adapter of a list in turn (fastp's `trimByMultiSequences`).
/// Returns the new length if at least one adapter was found.
pub fn trim_multi(seq: &[u8], adapters: &[Vec<u8>], p: &MatchParams) -> Option<usize> {
    let min_overlap = p.min_overlap_for(adapters.len());
    let mut len = seq.len();
    let mut trimmed = false;
    for adapter in adapters {
        if let Some(new_len) = find_adapter(&seq[..len], adapter, min_overlap, p) {
            len = new_len;
            trimmed = true;
        }
    }
    trimmed.then_some(len)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Test vectors from fastp's `AdapterTrimmer::test`.
    #[test]
    fn fastp_single_adapter() {
        let seq = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG";
        let len = trim_single(seq, b"TTTTCCACGGGGATACTACTG", &MatchParams::default()).unwrap();
        assert_eq!(&seq[..len], b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAA");
    }

    #[test]
    fn fastp_multi_adapter() {
        let seq = b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGGAAATTTCCCGGGAAATTTCCCGGGATCGATCGATCGATCGAATTCC";
        let adapters: Vec<Vec<u8>> = [
            "GCTAGCTAGCTAGCTA",
            "AAATTTCCCGGGAAATTTCCCGGG",
            "ATCGATCGATCGATCG",
            "AATTCCGGAATTCCGG",
        ]
        .iter()
        .map(|s| s.as_bytes().to_vec())
        .collect();
        let len = trim_multi(seq, &adapters, &MatchParams::default()).unwrap();
        assert_eq!(
            &seq[..len],
            b"TTTTAACCCCCCCCCCCCCCCCCCCCCCCCCCCCAATTTTAAAATTTTCCCCGGGG"
        );
    }

    #[test]
    fn adapter_dimer_empties_read() {
        // read starts with the adapter minus its first base (A-tailing)
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let seq = b"GATCGGAAGAGCACACGTCTGAACTCCAGTCAAAAA";
        assert_eq!(trim_single(seq, adapter, &MatchParams::default()), Some(0));
    }

    #[test]
    fn adapter_with_insertion_mid_read() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let insert = b"TTGCATGCAAGTCCATGGATCCATTGACCATTGGAC";
        let mut seq = insert.to_vec();
        seq.extend_from_slice(b"AGATCGGAAGTAGCACACGTCTGAACTCC"); // extra T after AGATCGGAAG
        assert_eq!(
            trim_single(&seq, adapter, &MatchParams::default()),
            Some(insert.len())
        );
    }

    #[test]
    fn no_adapter_found() {
        let seq = b"ACGTTGCAACGTTGCAACGTTGCAACGTTGCA";
        assert_eq!(
            trim_single(
                seq,
                b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
                &MatchParams::default()
            ),
            None
        );
    }

    #[test]
    fn mismatch_cap_is_respected() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let insert = b"TTGCATGCAAGTCCATGGATCCATTGACCATTGGAC";
        let mut seq = insert.to_vec();
        seq.extend_from_slice(b"AGATCGGAAGAGCACTCGTCTGAAC"); // one mismatch
        let strict = MatchParams {
            max_mismatches: Some(0),
            ..MatchParams::default()
        };
        assert_eq!(
            trim_single(&seq, adapter, &MatchParams::default()),
            Some(insert.len())
        );
        assert_ne!(trim_single(&seq, adapter, &strict), Some(insert.len()));
    }
}
