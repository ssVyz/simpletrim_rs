# Changelog

All notable changes to this project are documented here.
Format: [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), versioning: [SemVer](https://semver.org/).

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
