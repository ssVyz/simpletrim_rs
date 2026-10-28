//! Adapter auto-detection (port of fastp's `Evaluator::evalAdapterAndReadNum`).
//!
//! First the reads are scanned for the built-in adapter list. If no known adapter is frequent
//! enough, over-represented 10-mers are used as seeds and extended in both directions while
//! the reads agree, which assembles adapters that are not in the list.

use std::path::Path;

use crate::Result;
use crate::io::open_fastq;
use crate::known_adapters::KNOWN_ADAPTERS;
use crate::matcher::count_mismatches_bounded;

/// fastp reads up to 256K reads (and 151 bp per read on average) for detection.
pub const SAMPLE_READS: usize = 256 * 1024;
pub const SAMPLE_BASES: usize = 151 * SAMPLE_READS;
/// Fewer reads than this are not enough to detect anything.
pub const MIN_READS: usize = 10_000;
/// Detected adapters are cut to this length.
const MAX_ADAPTER_LEN: usize = 60;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Detected {
    pub sequence: String,
    /// Name from the built-in list, if the adapter is a known one.
    pub name: Option<&'static str>,
}

/// Outcome of adapter detection.
#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Detection {
    /// The adapter fastp would use.
    pub adapter: Option<Detected>,
    /// An over-represented sequence assembled from a seed k-mer that fastp's rules do not
    /// accept as an adapter (the extension did not end cleanly). Informational only.
    pub candidate: Option<String>,
}

/// Reads the sequences of the first reads of a file, within the fastp sampling limits.
pub fn load_sample(path: &Path, label: &str) -> Result<Vec<Vec<u8>>> {
    let mut input = open_fastq(path, label)?;
    let mut seqs = Vec::new();
    let mut bases = 0usize;
    while seqs.len() < SAMPLE_READS && bases < SAMPLE_BASES {
        let Some(chunk) = input.reader.read_chunk(4096, usize::MAX)? else {
            break;
        };
        for rec in chunk.iter() {
            if seqs.len() >= SAMPLE_READS || bases >= SAMPLE_BASES {
                break;
            }
            bases += rec.seq.len();
            seqs.push(rec.seq.to_vec());
        }
    }
    Ok(seqs)
}

/// Detects the adapter in a sample of read sequences. `trim_tail` is the configured fixed
/// tail trimming; the last cycle is always ignored because it is noisy.
pub fn detect(seqs: &[Vec<u8>], trim_tail: usize) -> Detection {
    if seqs.len() < MIN_READS {
        return Detection::default();
    }
    if let Some(known) = check_known_adapters(seqs) {
        return Detection {
            adapter: Some(known),
            candidate: None,
        };
    }
    detect_by_kmers(seqs, trim_tail.max(1))
}

/// Shorthand for the adapter of [`detect`].
pub fn detect_adapter(seqs: &[Vec<u8>], trim_tail: usize) -> Option<Detected> {
    detect(seqs, trim_tail).adapter
}

fn known(seq: &str, name: &'static str) -> Detected {
    Detected {
        sequence: seq.to_string(),
        name: Some(name),
    }
}

const KNOWN_MATCH_REQ: usize = 8;
const KNOWN_ONE_MISMATCH_PER: usize = 16;
/// Length of the exact seeds used to find candidate positions for known adapters.
const SEED_LEN: usize = 8;

/// fastp's test for a known adapter at `pos`: returns the mismatch count if it matches.
fn known_match_at(read: &[u8], adapter: &[u8], pos: usize) -> Option<usize> {
    let cmplen = (read.len() - pos).min(adapter.len());
    let allowed = cmplen / KNOWN_ONE_MISMATCH_PER;
    let mm = count_mismatches_bounded(&adapter[..cmplen], &read[pos..pos + cmplen], allowed);
    (mm <= allowed).then_some(mm)
}

fn seed_code(s: &[u8]) -> Option<u16> {
    s.iter()
        .try_fold(0u16, |k, &b| base_code(b).map(|c| (k << 2) | c as u16))
}

/// Calls `f(pos, code)` for every `SEED_LEN`-mer of `seq` without non-ACGT bases
/// (same encoding as [`seed_code`]).
fn for_each_seed(seq: &[u8], mut f: impl FnMut(usize, u16)) {
    let mut code = 0u16;
    let mut run = 0usize;
    for (i, &b) in seq.iter().enumerate() {
        match base_code(b) {
            Some(c) => {
                code = (code << 2) | c as u16;
                run += 1;
            }
            None => run = 0,
        }
        if run >= SEED_LEN {
            f(i + 1 - SEED_LEN, code);
        }
    }
}

/// Exact seeds of an adapter: `(offset, code)` of its first `len / 16 + 1` 8-mers.
///
/// A match at any position allows at most `cmplen / 16` mismatches over `cmplen >= 9` bases,
/// so at least one of the first `cmplen / 16 + 1` 8-mers of the adapter matches exactly
/// (pigeonhole). Checking only positions where a seed occurs therefore finds the same first
/// match as fastp's scan over every position.
fn adapter_seeds(adapter: &[u8]) -> Vec<(usize, u16)> {
    (0..=adapter.len() / KNOWN_ONE_MISMATCH_PER)
        .map(|j| j * SEED_LEN)
        .filter(|&off| off + SEED_LEN <= adapter.len())
        .filter_map(|off| seed_code(&adapter[off..off + SEED_LEN]).map(|c| (off, c)))
        .collect()
}

/// For every seed code, the `(adapter index, seed offset)` pairs of the known adapters.
struct SeedTable {
    start: Vec<u32>,
    entries: Vec<(u16, u8)>,
}

impl SeedTable {
    fn new() -> Self {
        let mut all: Vec<(u16, u16, u8)> = KNOWN_ADAPTERS
            .iter()
            .enumerate()
            .flat_map(|(i, (a, _))| {
                adapter_seeds(a.as_bytes())
                    .into_iter()
                    .map(move |(off, code)| (code, i as u16, off as u8))
            })
            .collect();
        all.sort_unstable();
        let mut start = vec![0u32; (1 << (2 * SEED_LEN)) + 1];
        for &(code, ..) in &all {
            start[code as usize + 1] += 1;
        }
        for i in 1..start.len() {
            start[i] += start[i - 1];
        }
        SeedTable {
            start,
            entries: all.into_iter().map(|(_, a, off)| (a, off)).collect(),
        }
    }

    fn get(&self, code: u16) -> &[(u16, u8)] {
        let c = code as usize;
        &self.entries[self.start[c] as usize..self.start[c + 1] as usize]
    }
}

/// Positions where a known adapter may match a read, as sorted `(adapter index, position)`
/// pairs: every position where one of the adapter's seeds occurs.
fn known_candidates(read: &[u8], table: &SeedTable, out: &mut Vec<(u16, usize)>) {
    out.clear();
    if read.len() <= KNOWN_MATCH_REQ {
        return;
    }
    let max_pos = read.len() - KNOWN_MATCH_REQ;
    for_each_seed(read, |p, code| {
        for &(adapter, off) in table.get(code) {
            let off = off as usize;
            if p >= off && p - off < max_pos {
                out.push((adapter, p - off));
            }
        }
    });
    out.sort_unstable();
    out.dedup();
}

/// Counts reads containing each known adapter (fastp `checkKnownAdapters`).
fn check_known_adapters(seqs: &[Vec<u8>]) -> Option<Detected> {
    const MAX_CHECK_READS: u64 = 100_000;
    const MAX_CHECK_BASES: u64 = MAX_CHECK_READS * 1000;
    const MAX_HIT: u64 = 1000;

    let table = SeedTable::new();
    let mut counts = vec![0u64; KNOWN_ADAPTERS.len()];
    let mut mismatches = vec![0u64; KNOWN_ADAPTERS.len()];
    let mut checked_reads = 0u64;
    let mut checked_bases = 0u64;
    let mut cur_max = 0u64;
    let mut candidates = Vec::new();

    for read in seqs {
        let rlen = read.len();
        checked_reads += 1;
        checked_bases += rlen as u64;
        if checked_reads > MAX_CHECK_READS || checked_bases > MAX_CHECK_BASES || cur_max > MAX_HIT {
            break;
        }
        known_candidates(read, &table, &mut candidates);
        // adapters in table order, as fastp iterates them; adapters without candidates
        // cannot match
        for group in candidates.chunk_by(|a, b| a.0 == b.0) {
            let i = group[0].0 as usize;
            let adapter = KNOWN_ADAPTERS[i].0.as_bytes();
            if adapter.len() >= rlen {
                continue;
            }
            // clearly not the winner, skip for speed
            if cur_max > 20 && counts[i] < cur_max / 10 {
                continue;
            }
            if let Some(mm) = group
                .iter()
                .find_map(|&(_, pos)| known_match_at(read, adapter, pos))
            {
                counts[i] += 1;
                cur_max = cur_max.max(counts[i]);
                mismatches[i] += mm as u64;
            }
        }
    }

    let mut best = None;
    let mut max_count = 0;
    for (i, &c) in counts.iter().enumerate() {
        if c > max_count {
            max_count = c;
            best = Some(i);
        }
    }
    let best = best?;
    let frequent = max_count > checked_reads / 50
        || (max_count > checked_reads / 200 && mismatches[best] < checked_reads);
    let (seq, name) = KNOWN_ADAPTERS[best];
    (frequent && seq.len() > 8).then(|| known(seq, name))
}

pub const KEYLEN: usize = 10;
pub const KEY_SPACE: usize = 1 << (KEYLEN * 2);

fn base_code(b: u8) -> Option<usize> {
    match b {
        b'A' => Some(0),
        b'T' | b'U' => Some(1),
        b'C' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// Rolling 2-bit encoding of `seq[pos..pos + KEYLEN]` (fastp `seq2int`); `last` is the key
/// at `pos - 1` or `None` to compute from scratch. Returns `None` if the window has an N.
pub fn seq2int(seq: &[u8], pos: usize, last: Option<usize>) -> Option<usize> {
    match last {
        Some(prev) => {
            let mask = KEY_SPACE - 1;
            base_code(seq[pos + KEYLEN - 1]).map(|c| ((prev << 2) & mask) + c)
        }
        None => seq[pos..pos + KEYLEN]
            .iter()
            .try_fold(0usize, |key, &b| base_code(b).map(|c| (key << 2) + c)),
    }
}

pub fn int2seq(mut key: usize) -> Vec<u8> {
    const BASES: [u8; 4] = [b'A', b'T', b'C', b'G'];
    let mut out = vec![b'N'; KEYLEN];
    for i in (0..KEYLEN).rev() {
        out[i] = BASES[key & 3];
        key >>= 2;
    }
    out
}

/// Calls `f(pos, key)` for every 10-mer starting at position 20 or later and ending at least
/// `shift_tail` bases before the read end (and starting before `max_pos`).
fn for_each_kmer(read: &[u8], shift_tail: usize, max_pos: usize, mut f: impl FnMut(usize, usize)) {
    let mut key = None;
    let mut pos = 20;
    while pos + KEYLEN + shift_tail <= read.len() && pos < max_pos {
        key = seq2int(read, pos, key);
        if let Some(k) = key {
            f(pos, k);
        }
        pos += 1;
    }
}

fn detect_by_kmers(seqs: &[Vec<u8>], shift_tail: usize) -> Detection {
    const TOP_N: usize = 10;
    const FOLD_THRESHOLD: u64 = 20;

    let mut counts = vec![0u32; KEY_SPACE];
    for read in seqs {
        for_each_kmer(read, shift_tail, usize::MAX, |_, k| counts[k] += 1);
    }
    // AAAAAAAAAA
    counts[0] = 0;

    // keep the TOP_N most frequent informative keys; ties go to the later key, as in fastp
    let mut top = [0usize; TOP_N];
    let mut total = 0u64;
    for k in 0..KEY_SPACE {
        let mut atcg = [0usize; 4];
        for i in 0..KEYLEN {
            atcg[(k >> (i * 2)) & 3] += 1;
        }
        let low_complexity = atcg.iter().any(|&n| n >= KEYLEN - 4);
        let too_much_gc = atcg[2] + atcg[3] >= KEYLEN - 2;
        let starts_with_gggg = k >> 12 == 0xff;
        if low_complexity || too_much_gc || starts_with_gggg {
            continue;
        }
        let val = counts[k];
        total += val as u64;
        for t in (0..TOP_N).rev() {
            if val < counts[top[t]] {
                if t < TOP_N - 1 {
                    top.copy_within(t + 1..TOP_N - 1, t + 2);
                    top[t + 1] = k;
                }
                break;
            } else if t == 0 {
                top.copy_within(0..TOP_N - 1, 1);
                top[0] = k;
            }
        }
    }

    // seeds worth extending, in rank order
    let mut seeds = Vec::new();
    for &key in &top {
        if key == 0 {
            continue;
        }
        let count = counts[key] as u64;
        if count < 10 || count * (KEY_SPACE as u64) < total * FOLD_THRESHOLD {
            break;
        }
        let seed = int2seq(key);
        let changes = seed.windows(2).filter(|w| w[0] != w[1]).count();
        if changes >= 3 {
            seeds.push(key);
        }
    }

    // (read, position) of every occurrence of each seed, collected in one pass
    const MAX_SEARCH_LENGTH: usize = 500;
    let mut seed_slot = vec![u8::MAX; KEY_SPACE];
    for (slot, &key) in seeds.iter().enumerate() {
        seed_slot[key] = slot as u8;
    }
    let mut occurrences: Vec<Vec<(usize, usize)>> = vec![Vec::new(); seeds.len()];
    if !seeds.is_empty() {
        for (r, read) in seqs.iter().enumerate() {
            for_each_kmer(read, shift_tail, MAX_SEARCH_LENGTH, |pos, k| {
                if let Some(occ) = occurrences.get_mut(seed_slot[k] as usize) {
                    occ.push((r, pos));
                }
            });
        }
    }

    let mut detection = Detection::default();
    for (&key, occ) in seeds.iter().zip(&occurrences) {
        match extend_seed(key, occ, seqs, shift_tail) {
            Ok(adapter) => {
                detection.adapter = Some(adapter);
                detection.candidate = None;
                return detection;
            }
            // a bare seed that could not be extended is not worth reporting
            Err(assembled) if assembled.len() >= 2 * KEYLEN => {
                detection.candidate.get_or_insert(assembled);
            }
            Err(_) => {}
        }
    }
    detection
}

/// Extends a seed 10-mer to a full adapter (fastp `getAdapterWithSeed`), given the
/// `(read, position)` occurrences of the seed. Returns the assembled sequence as `Err` if
/// fastp would reject it.
fn extend_seed(
    seed: usize,
    occurrences: &[(usize, usize)],
    seqs: &[Vec<u8>],
    shift_tail: usize,
) -> std::result::Result<Detected, String> {
    let all: Vec<usize> = (0..occurrences.len()).collect();
    let mut reached_leaf = true;
    // continuation after the seed, up to the noisy last cycles
    let forward_path = dominant_path(
        all.clone(),
        |i, depth| {
            let (r, pos) = occurrences[i];
            let read = &seqs[r];
            let idx = pos + KEYLEN + depth;
            (idx < read.len() - shift_tail).then(|| read[idx])
        },
        &mut reached_leaf,
    );
    // bases before the seed, walking backwards
    let backward_path = dominant_path(
        all,
        |i, depth| {
            let (r, pos) = occurrences[i];
            (depth < pos).then(|| seqs[r][pos - 1 - depth])
        },
        &mut reached_leaf,
    );

    let mut adapter: Vec<u8> = backward_path.into_iter().rev().collect();
    adapter.extend(int2seq(seed));
    adapter.extend(forward_path);
    adapter.truncate(MAX_ADAPTER_LEN);

    if let Some((seq, name)) = KNOWN_ADAPTERS
        .iter()
        .find(|(seq, _)| adapter.starts_with(seq.as_bytes()))
    {
        return Ok(known(seq, name));
    }
    let sequence = String::from_utf8_lossy(&adapter).into_owned();
    if reached_leaf {
        Ok(Detected {
            sequence,
            name: None,
        })
    } else {
        Err(sequence)
    }
}

/// Follows the dominant continuation of a set of sequences (fastp's
/// `NucleotideTree::getDominantPath`, computed without building the tree).
///
/// `base_at(i, depth)` gives base `depth` of sequence `i`, or `None` past its end; sequences
/// also stop at an `N`. At each depth the sequences still on the path vote with their next
/// base (bases are grouped by `base & 7` like fastp's tree slots); the walk continues while at
/// least 50 vote and one base gets 95% of the votes. Clears `reached_leaf` if it stopped
/// because no base was dominant.
fn dominant_path(
    mut active: Vec<usize>,
    base_at: impl Fn(usize, usize) -> Option<u8>,
    reached_leaf: &mut bool,
) -> Vec<u8> {
    const RATIO_THRESHOLD: f64 = 0.95;
    const NUM_THRESHOLD: u32 = 50;
    let slot_at = |i: usize, depth: usize| {
        base_at(i, depth)
            .filter(|&b| b != b'N')
            .map(|b| ((b & 0x07) as usize, b))
    };
    let mut path = Vec::new();
    loop {
        let depth = path.len();
        let mut counts = [0u32; 8];
        let mut first_base = [0u8; 8];
        for &i in &active {
            if let Some((slot, b)) = slot_at(i, depth) {
                if counts[slot] == 0 {
                    first_base[slot] = b;
                }
                counts[slot] += 1;
            }
        }
        let total: u32 = counts.iter().sum();
        if total < NUM_THRESHOLD {
            break;
        }
        let Some(slot) = (0..8).find(|&s| counts[s] as f64 / total as f64 >= RATIO_THRESHOLD)
        else {
            *reached_leaf = false;
            break;
        };
        path.push(first_base[slot]);
        active.retain(|&i| slot_at(i, depth).is_some_and(|(s, _)| s == slot));
    }
    path
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Deterministic pseudo-random insert sequences.
    fn random_seq(state: &mut u64, len: usize) -> Vec<u8> {
        (0..len)
            .map(|_| {
                *state = state
                    .wrapping_mul(6364136223846793005)
                    .wrapping_add(1442695040888963407);
                b"ACGT"[(*state >> 33) as usize % 4]
            })
            .collect()
    }

    fn reads_with_adapter(adapter: &[u8], n: usize) -> Vec<Vec<u8>> {
        let mut state = 42;
        (0..n)
            .map(|i| {
                let insert_len = if i % 3 == 0 { 60 + i % 50 } else { 150 };
                let mut read = random_seq(&mut state, insert_len.min(150));
                read.extend_from_slice(adapter);
                read.extend(random_seq(&mut state, 150));
                read.truncate(150);
                read
            })
            .collect()
    }

    #[test]
    fn dominant_path_matches_fastp_tree_test() {
        // same input as fastp's NucleotideTree::test
        let mut seqs: Vec<&[u8]> = Vec::new();
        for _ in 0..100 {
            seqs.extend([
                b"AAAATTTT".as_slice(),
                b"AAAATTTTGGGG",
                b"AAAATTTTGGGGCCCC",
                b"AAAATTTTGGGGCCAA",
            ]);
        }
        seqs.push(b"AAAATTTTGGGACCCC");
        let mut reached_leaf = true;
        let path = dominant_path(
            (0..seqs.len()).collect(),
            |i, d| seqs[i].get(d).copied(),
            &mut reached_leaf,
        );
        assert_eq!(path, b"AAAATTTTGGGGCC");
        assert!(!reached_leaf);
    }

    /// The seeded search must find exactly what fastp's scan over all positions finds.
    #[test]
    fn seeded_search_equals_full_scan() {
        let mut state = 1234;
        let table = SeedTable::new();
        let mut candidates = Vec::new();
        for round in 0..3000 {
            let planted_idx = round % KNOWN_ADAPTERS.len();
            let adapter = KNOWN_ADAPTERS[planted_idx].0.as_bytes();
            let insert = 20 + (round * 7) % 130;
            let mut read = random_seq(&mut state, insert);
            let mut planted = adapter.to_vec();
            // sprinkle mismatches into the planted adapter
            for k in 0..(round % 5) {
                let p = (k * 11 + round) % planted.len();
                planted[p] = b"ACGT"[(planted[p] as usize + 1) % 4];
            }
            read.extend(planted);
            read.extend(random_seq(&mut state, 150));
            read.truncate(150);
            if round % 13 == 0 {
                read[round % 150] = b'N';
            }
            known_candidates(&read, &table, &mut candidates);
            // every adapter, not only the planted one, must agree with the full scan
            for (i, (other, _)) in KNOWN_ADAPTERS.iter().enumerate() {
                let other = other.as_bytes();
                let fast = candidates
                    .iter()
                    .filter(|&&(a, _)| a as usize == i)
                    .find_map(|&(_, pos)| known_match_at(&read, other, pos).map(|mm| (pos, mm)));
                let full = (0..read.len() - KNOWN_MATCH_REQ)
                    .find_map(|p| known_match_at(&read, other, p).map(|mm| (p, mm)));
                assert_eq!(fast, full, "round {round}, adapter {i}");
                // at most one mismatch and at least 16 adapter bases in the read: must match
                if i == planted_idx && round % 5 < 2 && insert + 16 <= 150 {
                    assert!(full.is_some(), "planted adapter not found in round {round}");
                }
            }
        }
    }

    #[test]
    fn kmer_encoding_roundtrip() {
        let s = b"ATCGATCGAT";
        assert_eq!(int2seq(seq2int(s, 0, None).unwrap()), s);
        assert_eq!(seq2int(b"ATCGNTCGAT", 0, None), None);
    }

    #[test]
    fn detects_known_adapter() {
        let adapter = b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA";
        let reads = reads_with_adapter(adapter, 12_000);
        let d = detect_adapter(&reads, 0).unwrap();
        assert_eq!(d.sequence.as_bytes(), adapter);
        assert_eq!(d.name, Some("Illumina TruSeq Adapter Read 1"));
    }

    #[test]
    fn unknown_adapter_followed_by_random_bases_is_only_a_candidate() {
        let adapter = b"TCAGTGCATCAGTACGATCAGACTGACATGCAGTCAG";
        let reads = reads_with_adapter(adapter, 12_000);
        let d = detect(&reads, 0);
        assert_eq!(d.adapter, None);
        assert_eq!(
            d.candidate.as_deref(),
            Some(std::str::from_utf8(adapter).unwrap())
        );
    }

    #[test]
    fn nothing_in_clean_reads() {
        let mut state = 7;
        let reads: Vec<Vec<u8>> = (0..12_000).map(|_| random_seq(&mut state, 150)).collect();
        assert_eq!(detect_adapter(&reads, 0), None);
    }

    #[test]
    fn too_few_reads() {
        let reads = reads_with_adapter(b"AGATCGGAAGAGCACACGTCTGAACTCCAGTCA", 100);
        assert_eq!(detect_adapter(&reads, 0), None);
    }
}
