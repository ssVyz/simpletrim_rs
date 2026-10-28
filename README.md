# simpletrim

A small, portable NGS read trimmer written in Rust. It is configured entirely from the command
line (no config files, no environment variables), which makes it easy to call from other
programs.

**Status:** stage 2. Short reads (Illumina-style, single-end and paired-end) and long reads
(Oxford Nanopore, PacBio) are supported. Long-read trimming has so far been tested on synthetic
reads and one small nanopore run.

## AI disclosure

This program was created using agentic coding via claude code.

## Build

```sh
cargo build --release   # binary: target/release/simpletrim(.exe)
```

## Usage

```sh
# Paired-end trimming (.gz in/out is chosen by file extension)
simpletrim short --r1 in_R1.fastq.gz --r2 in_R2.fastq.gz --o1 out_R1.fastq.gz --o2 out_R2.fastq.gz

# Single-end trimming
simpletrim short --r1 in.fastq.gz --o1 out.fastq.gz

# Long reads (nanopore, PacBio)
simpletrim long -i reads.fastq.gz -o trimmed.fastq.gz

# Inspect the input and get a suggested `simpletrim short` or `simpletrim long` command
simpletrim analyze --r1 in_R1.fastq.gz --r2 in_R2.fastq.gz
simpletrim analyze --r1 reads.fastq.gz
```

`simpletrim short` applies these steps in order, and each one can be switched off:

1. Fixed front/tail trimming.
2. Sliding-window quality cutting (`--cut-front`, `--cut-tail`, `--cut-right`; off by default).
3. Poly-G tail trimming (on automatically for two-colour instruments).
4. Adapter trimming: read overlap for pairs, then adapter sequences.
   - Adapters can be built-in known adapters, presets, custom sequences, a FASTA list, or
     auto-detected.
5. Adapter-dimer check (drops the whole pair).
6. Quality and length filters.

For paired-end data, mate names are checked. `--orphan-mode discard|keep-both|orphan` decides
what happens when only one mate passes. Run `simpletrim short --help` for all options and
defaults.

`simpletrim long` applies these steps in order:

1. Fixed front/tail trimming.
2. Sliding-window quality cutting (`--cut-front`, `--cut-tail`; off by default).
3. Poly-X tail trimming (`--trim-poly-x`; off by default).
4. Adapter trimming at the read start and end, also when the read starts or ends inside the
   adapter. Adapters are auto-detected unless given (`--start-adapter`, `--end-adapter`,
   `--adapter-fasta`). Matching uses edit distance, so insertions and deletions are tolerated.
5. Chimeric reads: reads are split where an adapter is found inside them
   (`--discard-chimeric-reads` drops them instead, `--disable-chimera-splitting` keeps them).
6. Low-quality stretches inside reads (`--low-quality-mode split|longest`; off by default):
   `split` writes every good piece as a read of its own, `longest` keeps only the longest
   piece.
7. Quality, length and optional GC and complexity filters, applied to every piece.

Pieces of a split read are named `<id>_1`, `<id>_2`, ... (read order); the rest of the header
is kept. `--failed-output FILE` collects reads and pieces that fail, with the reason appended
to the name. Quality means average error probabilities, as is usual for nanopore data.

### Calling from other programs

- `--json` writes progress and results as one JSON object per line on stderr. The `analyze`
  report is printed on stdout (as JSON with `--json`).
- `--report FILE` saves the run summary of `short` or `long` as JSON.
- Exit codes: `0` success, `1` error, `2` invalid command line.
- Existing output files are never overwritten unless `--force` is given. Outputs of a failed
  run are removed.
- `-w/--threads` sets the number of worker threads (default 3).

### Synthetic test data

The long-read tests use synthetic nanopore-like reads with a known structure
(`tests/common/synth.rs`). A larger reference set (about 185 Mbp, with a truth table of insert
and low-quality positions) can be written to `examples/` with:

```sh
cargo test --release --test long reference_dataset -- --ignored
```

## Acknowledgements

### fastp

simpletrim was written using [fastp](https://github.com/OpenGene/fastp) (version 1.4.0) as a
template. Its algorithms, defaults and command-line behaviour were taken from fastp's source
code and reimplemented in Rust. This includes:

- fixed and sliding-window trimming;
- poly-G trimming;
- adapter detection and trimming by sequence and by read overlap;
- adapter-dimer handling;
- the quality and length filters.

The built-in adapter list in `src/known_adapters.rs` is taken from fastp's `knownadapters.h`,
which in turn draws on [stephenturner/adapters](https://github.com/stephenturner/adapters).
The deliberate differences from fastp are listed in [CHANGELOG.md](CHANGELOG.md).

simpletrim is not affiliated with or endorsed by the fastp authors.

Many thanks to Shifu Chen, the OpenGene contributors and HaploX Biotechnology for developing
and maintaining fastp. If simpletrim is useful in your work, please also cite fastp:

- Shifu Chen. fastp 1.0: An ultra-fast all-round tool for FASTQ data quality control and
  preprocessing. *iMeta* 4.5 (2025): e70078. <https://doi.org/10.1002/imt2.70078>
- Shifu Chen. Ultrafast one-pass FASTQ data preprocessing, quality control, and deduplication
  using fastp. *iMeta* 2 (2023): e107. <https://doi.org/10.1002/imt2.107>
- Shifu Chen, Yanqing Zhou, Yaru Chen, Jia Gu. fastp: an ultra-fast all-in-one FASTQ
  preprocessor. *Bioinformatics* 34.17 (2018): i884–i890.
  <https://doi.org/10.1093/bioinformatics/bty560>

### fastplong

`simpletrim long` was written using [fastplong](https://github.com/OpenGene/fastplong)
(version 0.7.1, MIT License, by Shifu Chen and the OpenGene contributors) as a template, in
the same way: its processing steps, adapter auto-detection, chimera splitting, low-quality
segment handling, filters and defaults were taken from its source code and reimplemented in
Rust. simpletrim is not affiliated with or endorsed by the fastplong authors. The deliberate
differences (among them edit-distance adapter matching and quality means from error
probabilities) are listed in [CHANGELOG.md](CHANGELOG.md).

### Edit distance

Adapters in long reads are found with the bit-parallel edit-distance algorithm of Gene Myers
(in Heikki Hyyrö's formulation), implemented from the publication:

- Gene Myers. A fast bit-vector algorithm for approximate string matching based on dynamic
  programming. *Journal of the ACM* 46.3 (1999): 395–415. <https://doi.org/10.1145/316542.316550>

### Rust crates

simpletrim depends on these crates. Thanks to their authors and maintainers.

| Crate | Used for | License |
|---|---|---|
| [clap](https://github.com/clap-rs/clap) | command-line parsing | MIT OR Apache-2.0 |
| [flate2](https://github.com/rust-lang/flate2-rs) | gzip reading and writing | MIT OR Apache-2.0 |
| [zlib-rs](https://github.com/trifectatechfoundation/zlib-rs) | pure-Rust zlib backend of flate2 | Zlib |
| [rayon](https://github.com/rayon-rs/rayon) | parallel processing | MIT OR Apache-2.0 |
| [serde](https://github.com/serde-rs/serde) | serialization | MIT OR Apache-2.0 |
| [serde_json](https://github.com/serde-rs/json) | JSON output | MIT OR Apache-2.0 |

Their own dependencies are listed in `Cargo.lock` (`cargo tree` shows the full tree).
