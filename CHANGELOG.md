# Changelog

All notable changes to this project are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioning: [SemVer](https://semver.org/).

## [0.0.3] - 2026-10-28

### Added
- `simpletrim long`: trimming of long reads (Oxford Nanopore, PacBio), following fastplong's
  steps and defaults. Steps, in order: fixed front/tail trimming, sliding-window quality
  cutting (`--cut-front/--cut-tail`, off by default), poly-X tail trimming (`--trim-poly-x`,
  off by default), adapter trimming at the read start and end (`--start-adapter`,
  `--end-adapter`, `--adapter-fasta`; auto-detected by default), splitting of chimeric reads
  at internal adapters (`--discard-chimeric-reads`, `--disable-chimera-splitting`), removal of
  low-quality stretches (`--low-quality-mode split|longest`, off by default), and quality,
  N, length, GC and complexity filters applied to every piece.
  - Pieces of a split read are written as `<id>_1`, `<id>_2`, ... with the rest of the header
    kept; `--failed-output` collects failing reads and pieces with fastplong's reason tags.
  - JSON events and `--report` as for `short`, with read length (N50, N90), read quality and
    per-step statistics.
- `simpletrim analyze` for long reads (detected from the header, or a median read length of
  at least 1000 bp): lengths and N50, read qualities, quality along the reads, nanopore run
  details from Dorado/MinKNOW headers, adapter detection and how often adapters occur at the
  ends and inside reads, poly-X tails and low-quality stretches; suggests a `simpletrim long`
  command. New option `--max-bases` (default 500 Mbp).
- `simpletrim short`: poly-X tail trimming after adapter trimming (`--trim-poly-x`,
  `--poly-x-min-len`, as fastp).
- Synthetic nanopore-like test reads (`tests/common/synth.rs`) and an ignored test that writes a
  reference dataset with a truth table to `examples/`.

### Changed
- The reader/worker/writer pipeline is shared by `short` and `long`; chunks of long reads are
  limited by size (8 MB) rather than only by record count.

### Differences from fastplong
- Adapters are matched by edit distance (Myers' bit-parallel algorithm), so insertions and
  deletions are tolerated anywhere; fastplong compares position by position and checks the
  edit distance of only one candidate position.
- All adapters of one read end are searched in the same window of the untrimmed read and the
  innermost match decides the cut; fastplong applies them one after another and moves its
  window after every cut.
- Adapters inside reads use a stricter error rate (`--chimera-distance-threshold`, default
  0.15) and must be at least 20 bases long: on random sequence an error rate of 0.25 gives
  about 8 false matches per Mbp for a 24 bp adapter, 0.15 about 0.25 (calibration test in
  `src/edit.rs`).
- Quality means average error probabilities (sliding windows, low-quality stretches, the
  `--average-quality` filter and all statistics). Default thresholds were converted so that
  they flag as many windows on real R10.4.1 data as fastplong's plain-mean thresholds:
  `--cut-mean-quality` 17 (fastplong 20), `--low-quality-threshold` 7 (fastplong 10).
- Low-quality stretches are found correctly in the whole read; fastplong's detection (used by
  `--break` and `--mask`) loses its running sum after the first stretch and marks everything
  from there to the read end. Stretch edges are shrunk to the first and last poor base.
- `--low-quality-mode longest` keeps the longest good piece; fastplong's `--best_read_segment`
  keeps the piece with the best mean quality, which can be a very short one. Masking with N
  (`--mask`) is not offered.
- Adapter detection skips 10-mers with fewer than three base changes, as intended; fastplong
  applies this test to the 10-mer's count by mistake.
- Pieces of split reads keep their read id with a suffix (`<id>_1`); fastplong prefixes the
  name (`chimera-segment-1-<id>`, `r1-<id>`).
- A read that is trimmed away completely is counted as too short and can be written to
  `--failed-output`; fastplong drops it without counting it.
- Existing outputs are never overwritten without `--force` (fastplong overwrites by default).
- Not included: sampling (`--sample_rate`, `--target_bases`), best-read selection
  (`--best_reads`, `--best_bases`), output splitting, `--stdin/--stdout`, ONT sequencing
  summary parsing, HTML reports.

## [0.0.2] - 2026-10-28

### Added
- `simpletrim short`: trimming of Illumina-style single-end and paired-end reads (plain or
  `.gz` FASTQ in and out), following fastp's algorithms and defaults. Steps, in order, each
  can be switched off: fixed front/tail trimming, sliding-window quality cutting
  (`--cut-front/--cut-tail/--cut-right`, off by default), poly-G trimming (automatic for
  two-colour instruments), adapter trimming (read overlap for pairs, then adapter sequences),
  adapter dimer check, quality filter (on by default, as fastp) and length filter.
- Adapters: built-in list of 234 known adapters (from fastp), presets (`--adapter-preset
  truseq|nextera|smallrna|mgi`), custom sequences, FASTA lists and fastp-style
  auto-detection (single-end by default, paired-end with `--detect-adapter-pe`). Overlap and
  sequence matching tolerances (minimum overlap, mismatches, error rate) are configurable.
- Paired-end: read names are checked (mismatch aborts unless `--no-name-check`), unequal
  read counts are an error, and `--orphan-mode discard|keep-both|orphan` decides what happens
  when only one mate passes. Adapter dimers always drop the whole pair.
- `simpletrim analyze`: inspects a sample of the input (default 1M reads/pairs) and reports
  platform/instrument from headers, pairing, lengths, qualities, GC, N, poly-G, adapter
  motifs, composition bias, read overlap/insert size and detected adapters, then suggests a
  `simpletrim short` command. Long-read (ONT/PacBio) input is recognised; its analysis is a
  stub.
- `simpletrim long`: placeholder, reports "not implemented yet".
- Human readable progress on stderr, or JSON lines with `--json`; `--quiet`; `--report` writes
  the run summary as JSON. Exit codes: 0 success, 1 error, 2 invalid command line.
- Streaming chunked processing on a rayon worker pool (`--threads`, default 3); `.gz` output
  is compressed in parallel as multi-member gzip (`--compression-level`, default 4).
- Existing outputs are never overwritten without `--force`; outputs of a failed run are
  removed.
- Dependencies: clap, flate2 (pure-Rust zlib-rs backend), rayon, serde, serde_json.
- `README.md` with usage, acknowledgement of the Rust crates used, and disclosure of fastp as
  the template for this program (authors, citations, MIT license notice).

### Differences from fastp
- Bug fixes: gapped (one insertion/deletion) adapter matching compares at the candidate
  position instead of the read start; gapped overlap analysis no longer rejects valid
  alignments; fixed trimming combined with `--cut-front/--cut-tail` no longer drops
  `window - 1` extra good bases; `--cut-right` never keeps fixed-trimmed tail bases.
- Poly-G auto-detection requires a digit after `@A` (NovaSeq 6000), so Element AVITI (`@AV`)
  reads are not treated as two-colour data.
- When one mate is fully trimmed away, the other mate still gets poly-G and adapter trimming
  (it matters for `--orphan-mode`).
- Not included: merging, deduplication, UMI handling, base correction, poly-X, low complexity
  and index filters, phred64 input, HTML reports.

## [0.0.1] - 2026-10-27

### Added
- Initial project setup: Cargo binary crate, `.gitignore`, `CHANGELOG.md`, `CLAUDE.md`.
- Binary is named `simpletrim` (builds `simpletrim.exe`); package name stays `simpletrim_rs`.
