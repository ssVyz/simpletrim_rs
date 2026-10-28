//! Approximate matching of adapters in reads with edit distance (substitutions, insertions
//! and deletions), as needed for nanopore reads where insertions and deletions are the main
//! error type.
//!
//! Patterns of up to 64 bases use Myers' bit-parallel algorithm (G. Myers, J. ACM 46(3),
//! 1999, in Hyyrö's formulation): one pass over the text with a few integer operations per
//! base, no SIMD. Longer patterns fall back to the textbook dynamic programme, which is also
//! the reference in the tests.
//!
//! Bases are compared case-insensitively, `U` equals `T`, and an `N` in the pattern matches
//! any base. An `N` (or any other character) in the text only matches an `N` in the pattern.

/// Base code used for matching: A, C, G, T/U, and 4 for everything else.
const fn code_table() -> [u8; 256] {
    let mut t = [4u8; 256];
    t[b'A' as usize] = 0;
    t[b'a' as usize] = 0;
    t[b'C' as usize] = 1;
    t[b'c' as usize] = 1;
    t[b'G' as usize] = 2;
    t[b'g' as usize] = 2;
    t[b'T' as usize] = 3;
    t[b't' as usize] = 3;
    t[b'U' as usize] = 3;
    t[b'u' as usize] = 3;
    t
}

static CODE: [u8; 256] = code_table();

/// A pattern prepared for searching.
pub struct Pattern {
    seq: Vec<u8>,
    /// Bit `i` of `peq[c]` is set if pattern base `i` matches a text base with code `c`.
    peq: [u64; 5],
}

impl Pattern {
    pub fn new(seq: &[u8]) -> Self {
        let mut peq = [0u64; 5];
        if seq.len() <= 64 {
            for (i, &b) in seq.iter().enumerate() {
                let code = CODE[b as usize];
                if code == 4 {
                    // N matches anything
                    for p in &mut peq {
                        *p |= 1 << i;
                    }
                } else {
                    peq[code as usize] |= 1 << i;
                }
            }
        }
        Pattern {
            seq: seq.to_vec(),
            peq,
        }
    }

    pub fn len(&self) -> usize {
        self.seq.len()
    }

    pub fn seq(&self) -> &[u8] {
        &self.seq
    }
}

/// Whether pattern base `p` matches text base `t`.
fn matches(p: u8, t: u8) -> bool {
    let (cp, ct) = (CODE[p as usize], CODE[t as usize]);
    cp == 4 || (cp == ct && ct != 4)
}

/// Calls `f(j, d)` for every text position `j`, where `d` is the smallest edit distance of the
/// whole pattern to a substring of `text` that ends at `j` (inclusive).
///
/// With `overhang`, the pattern may also start before the text: a prefix of the pattern can
/// hang off the start of the text without cost. Then `d` at position `j` also covers
/// alignments of a pattern suffix to `text[..=j]`.
pub fn for_each_end(p: &Pattern, text: &[u8], overhang: bool, f: impl FnMut(usize, u32)) {
    let m = p.len();
    if m == 0 {
        return;
    }
    if m <= 64 {
        myers(p, text, overhang, f);
    } else {
        dp(p.seq(), text, overhang, f);
    }
}

fn myers(p: &Pattern, text: &[u8], overhang: bool, mut f: impl FnMut(usize, u32)) {
    let m = p.len();
    let high = 1u64 << (m - 1);
    // vertical deltas of the first column: +1 per row (D[i][0] = i), or 0 with overhang
    let (mut pv, mut mv, mut score) = if overhang {
        (0u64, 0u64, 0u32)
    } else {
        (u64::MAX, 0u64, m as u32)
    };
    for (j, &c) in text.iter().enumerate() {
        let eq = p.peq[CODE[c as usize] as usize];
        let xv = eq | mv;
        let xh = ((eq & pv).wrapping_add(pv) ^ pv) | eq;
        let mut ph = mv | !(xh | pv);
        let mut mh = pv & xh;
        if ph & high != 0 {
            score += 1;
        } else if mh & high != 0 {
            score -= 1;
        }
        // the first row is 0 everywhere (the match may start anywhere in the text)
        ph <<= 1;
        mh <<= 1;
        pv = mh | !(xv | ph);
        mv = ph & xv;
        f(j, score);
    }
}

/// Dynamic programme with the same semantics as [`myers`], for any pattern length.
fn dp(pattern: &[u8], text: &[u8], overhang: bool, mut f: impl FnMut(usize, u32)) {
    let m = pattern.len();
    let mut col: Vec<u32> = (0..=m as u32)
        .map(|i| if overhang { 0 } else { i })
        .collect();
    for (j, &c) in text.iter().enumerate() {
        let mut diag = col[0];
        col[0] = 0;
        for i in 1..=m {
            let up = col[i];
            let sub = diag + !matches(pattern[i - 1], c) as u32;
            col[i] = sub.min(up + 1).min(col[i - 1] + 1);
            diag = up;
        }
        f(j, col[m]);
    }
}

/// A match: the text position of its last base and its edit distance.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Hit {
    pub end: usize,
    pub dist: u32,
}

/// Allowed edit distance for an aligned length (`round(rate * len)`, as fastplong).
pub fn allowed(rate: f64, len: usize) -> u32 {
    (rate * len as f64).round() as u32
}

/// Finds the pattern at the start of `text`, where it may be cut off by the text start: either
/// the whole pattern anywhere in `text`, or a suffix of at least `min_overlap` pattern bases
/// that starts at the text start. The allowed distance is `rate` times the aligned length
/// (the pattern length, or the number of text bases for a cut-off pattern).
///
/// Of the matching positions, the innermost (rightmost) group of consecutive positions is
/// taken, and in it the position with the smallest distance (the rightmost on ties).
pub fn find_innermost(p: &Pattern, text: &[u8], rate: f64, min_overlap: usize) -> Option<Hit> {
    let m = p.len();
    let mut best: Option<Hit> = None;
    let mut last_hit_end = None;
    for_each_end(p, text, true, |j, d| {
        let aligned = (j + 1).min(m);
        if j + 1 < min_overlap || d > allowed(rate, aligned) {
            return;
        }
        let continues = last_hit_end == Some(j.wrapping_sub(1));
        last_hit_end = Some(j);
        match best {
            Some(b) if continues && b.dist < d => {}
            _ => best = Some(Hit { end: j, dist: d }),
        }
    });
    best
}

/// Finds all occurrences of the whole pattern with at most `max_dist` edits. Consecutive
/// matching positions form one occurrence, reported at its smallest distance (the leftmost on
/// ties).
pub fn find_all(p: &Pattern, text: &[u8], max_dist: u32) -> Vec<Hit> {
    let mut hits: Vec<Hit> = Vec::new();
    let mut last_hit_end = None;
    for_each_end(p, text, false, |j, d| {
        if d > max_dist {
            return;
        }
        let continues = last_hit_end == Some(j.wrapping_sub(1));
        last_hit_end = Some(j);
        match hits.last_mut() {
            Some(h) if continues => {
                if d < h.dist {
                    *h = Hit { end: j, dist: d };
                }
            }
            _ => hits.push(Hit { end: j, dist: d }),
        }
    });
    hits
}

#[cfg(test)]
mod tests {
    use super::*;

    struct Lcg(u64);

    impl Lcg {
        fn next(&mut self) -> u64 {
            self.0 = self
                .0
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            self.0 >> 33
        }

        fn below(&mut self, n: usize) -> usize {
            (self.next() % n as u64) as usize
        }

        fn seq(&mut self, len: usize, alphabet: &[u8]) -> Vec<u8> {
            (0..len)
                .map(|_| alphabet[self.below(alphabet.len())])
                .collect()
        }
    }

    fn ends(p: &[u8], text: &[u8], overhang: bool, use_dp: bool) -> Vec<u32> {
        let mut out = Vec::new();
        if use_dp {
            dp(p, text, overhang, |_, d| out.push(d));
        } else {
            for_each_end(&Pattern::new(p), text, overhang, |_, d| out.push(d));
        }
        out
    }

    /// Plain edit distance of two strings, for checking the semantics.
    fn edit_distance(a: &[u8], b: &[u8]) -> u32 {
        let mut prev: Vec<u32> = (0..=b.len() as u32).collect();
        for (i, &x) in a.iter().enumerate() {
            let mut cur = vec![i as u32 + 1; b.len() + 1];
            for (j, &y) in b.iter().enumerate() {
                cur[j + 1] = (prev[j] + !matches(x, y) as u32)
                    .min(prev[j + 1] + 1)
                    .min(cur[j] + 1);
            }
            prev = cur;
        }
        prev[b.len()]
    }

    #[test]
    fn myers_equals_dp() {
        let mut rng = Lcg(7);
        for round in 0..3000 {
            let m = 1 + rng.below(64);
            let n = rng.below(150);
            // small alphabets give many near matches; N and lowercase exercise the rules
            let alphabet: &[u8] = match round % 4 {
                0 => b"AC",
                1 => b"ACGT",
                2 => b"ACGTNacgu",
                _ => b"AAAC",
            };
            let p = rng.seq(m, alphabet);
            let t = rng.seq(n, alphabet);
            for overhang in [false, true] {
                assert_eq!(
                    ends(&p, &t, overhang, false),
                    ends(&p, &t, overhang, true),
                    "pattern {} text {} overhang {overhang}",
                    String::from_utf8_lossy(&p),
                    String::from_utf8_lossy(&t)
                );
            }
        }
    }

    #[test]
    fn semantics_match_plain_edit_distance() {
        let mut rng = Lcg(11);
        for _ in 0..300 {
            let (m, n) = (1 + rng.below(12), rng.below(14));
            let p = rng.seq(m, b"ACG");
            let t = rng.seq(n, b"ACG");
            let free = ends(&p, &t, false, false);
            let hang = ends(&p, &t, true, false);
            for j in 0..t.len() {
                // whole pattern against any substring ending at j
                let best = (0..=j + 1)
                    .map(|s| edit_distance(&p, &t[s..=j]))
                    .min()
                    .unwrap();
                assert_eq!(free[j], best);
                // ... or a pattern suffix against a text prefix
                let suffix = (0..=p.len())
                    .map(|k| edit_distance(&p[k..], &t[..=j]))
                    .min()
                    .unwrap();
                assert_eq!(hang[j], best.min(suffix));
            }
        }
    }

    #[test]
    fn long_patterns_use_dp() {
        let mut rng = Lcg(3);
        let p = rng.seq(80, b"ACGT");
        let mut t = rng.seq(50, b"ACGT");
        t.extend_from_slice(&p);
        t.extend(rng.seq(50, b"ACGT"));
        let hits = find_all(&Pattern::new(&p), &t, 5);
        assert_eq!(hits, vec![Hit { end: 129, dist: 0 }]);
    }

    #[test]
    fn rules_for_n_u_and_case() {
        let p = Pattern::new(b"ACGNT");
        assert_eq!(find_all(&p, b"xxacgaUyy", 0), vec![Hit { end: 6, dist: 0 }]);
        // N in the text does not match
        assert_eq!(find_all(&Pattern::new(b"ACGT"), b"ACNT", 0), vec![]);
    }

    #[test]
    fn innermost_partial_and_full() {
        let adapter = b"GTTCAGTTACGTATTGCTAAGGTTAA";
        let p = Pattern::new(adapter);
        // read starts with the last 18 adapter bases
        let mut read = adapter[8..].to_vec();
        read.extend_from_slice(b"CCCCATTTGGGACCCATTTGACCCA");
        let h = find_innermost(&p, &read, 0.25, 16).unwrap();
        assert_eq!((h.end, h.dist), (17, 0));
        // with one deletion inside the adapter part
        let mut del = read.clone();
        del.remove(5);
        let h = find_innermost(&p, &del, 0.25, 16).unwrap();
        assert_eq!((h.end, h.dist), (16, 1));
        // too little overlap
        assert_eq!(find_innermost(&p, &read[10..], 0.25, 16), None);
        // two full copies: the inner one wins
        let mut two = b"TTT".to_vec();
        two.extend_from_slice(adapter);
        two.extend_from_slice(b"GGGGG");
        two.extend_from_slice(adapter);
        two.extend_from_slice(b"CACACACACACA");
        let h = find_innermost(&p, &two, 0.25, 16).unwrap();
        assert_eq!((h.end, h.dist), (3 + 26 + 5 + 26 - 1, 0));
    }

    #[test]
    fn find_all_reports_each_occurrence_once() {
        let a = b"ACGTTGCAAGGCTTAACCGGTTAA";
        let mut t = b"CCCCCCCCCC".to_vec();
        t.extend_from_slice(a);
        t.extend_from_slice(b"TTTTTTTTTTTTTTTTTTTT");
        let mut b = a.to_vec();
        b.insert(10, b'G'); // one insertion
        t.extend_from_slice(&b);
        t.extend_from_slice(b"CCCCC");
        let hits = find_all(&Pattern::new(a), &t, 3);
        assert_eq!(hits.len(), 2);
        assert_eq!(hits[0], Hit { end: 33, dist: 0 });
        assert_eq!(hits[1], Hit { end: 78, dist: 1 });
    }

    /// Prints false-positive and detection rates of the matching rules on synthetic data.
    /// Run with `cargo test --release calibration -- --ignored --nocapture`.
    #[test]
    #[ignore]
    fn calibration_report() {
        use crate::synth::{Errors, Rng, mutate_into};
        let mut rng = Rng::new(42);
        let text = rng.seq(4_000_000);
        println!("random hits of the whole pattern (find_all) per Mbp of random sequence");
        for m in [20usize, 24, 30, 40, 50, 64] {
            let p = Pattern::new(&rng.seq(m));
            let mut row = format!("  m={m:2}:");
            for rate in [0.10, 0.15, 0.20, 0.25] {
                let hits = find_all(&p, &text, allowed(rate, m)).len();
                row.push_str(&format!("  rate {rate:.2}: {:9.3}", hits as f64 / 4.0));
            }
            println!("{row}");
        }
        println!("detection of a whole adapter with nanopore-like errors inside random sequence");
        for m in [24usize, 40, 64] {
            let adapter = rng.seq(m);
            let p = Pattern::new(&adapter);
            for err in [0.05, 0.10, 0.15] {
                let mut row = format!("  m={m:2} errors {err:.2}:");
                for rate in [0.10, 0.15, 0.20, 0.25] {
                    let mut found = 0;
                    for _ in 0..2000 {
                        let mut seq = rng.seq(100);
                        let mut qual = Vec::new();
                        mutate_into(
                            &mut rng,
                            &adapter,
                            Errors::nanopore(err),
                            (20, 20),
                            &mut seq,
                            &mut qual,
                        );
                        seq.extend(rng.seq(100));
                        if !find_all(&p, &seq, allowed(rate, m)).is_empty() {
                            found += 1;
                        }
                    }
                    row.push_str(&format!("  rate {rate:.2}: {:5.1}%", found as f64 / 20.0));
                }
                println!("{row}");
            }
        }
        println!("read starts: random 200 bp windows with a hit (find_innermost, rate 0.25)");
        for m in [24usize, 40, 64] {
            let p = Pattern::new(&rng.seq(m));
            let mut row = format!("  m={m:2}:");
            for min_overlap in [10usize, 12, 14, 16, 20] {
                let mut hits = 0;
                for _ in 0..20000 {
                    if find_innermost(&p, &rng.seq(200), 0.25, min_overlap).is_some() {
                        hits += 1;
                    }
                }
                row.push_str(&format!(
                    "  min overlap {min_overlap}: {:6.3}%",
                    hits as f64 / 200.0
                ));
            }
            println!("{row}");
        }
        println!(
            "read starts: detection of a cut-off adapter (last L bases) with 6% / 12% errors, rate 0.25, min overlap 16"
        );
        let adapter = rng.seq(40);
        let p = Pattern::new(&adapter);
        for l in [12usize, 16, 20, 25, 30, 40] {
            let mut row = format!("  L={l:2}:");
            for err in [0.06, 0.12] {
                let (mut found, mut exact) = (0, 0);
                for _ in 0..2000 {
                    let (mut seq, mut qual) = (Vec::new(), Vec::new());
                    mutate_into(
                        &mut rng,
                        &adapter[40 - l..],
                        Errors::nanopore(err),
                        (20, 20),
                        &mut seq,
                        &mut qual,
                    );
                    let true_end = seq.len();
                    seq.extend(rng.seq(200));
                    if let Some(h) = find_innermost(&p, &seq[..200], 0.25, 16) {
                        found += 1;
                        if (h.end + 1).abs_diff(true_end) <= 2 {
                            exact += 1;
                        }
                    }
                }
                row.push_str(&format!(
                    "  errors {err:.2}: found {:5.1}% (end within 2 bp: {:5.1}%)",
                    found as f64 / 20.0,
                    exact as f64 / 20.0
                ));
            }
            println!("{row}");
        }
    }
}
