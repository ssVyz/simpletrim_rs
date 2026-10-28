//! Command line interface.

use std::path::PathBuf;

use clap::{Args, Parser, Subcommand, ValueEnum};

use crate::known_adapters::ADAPTER_PRESETS;

#[derive(Parser, Debug)]
#[command(
    name = "simpletrim",
    bin_name = "simpletrim",
    version,
    about = "A small, portable NGS read trimmer",
    long_about = "A small, portable NGS read trimmer.\n\n\
        All settings are given as command line arguments; there are no config files or \
        environment variables. Progress and results are reported on stderr, as text or with \
        --json as one JSON object per line.\n\n\
        Exit codes: 0 success, 1 error, 2 invalid command line.",
    propagate_version = true
)]
pub struct Cli {
    /// Machine-readable output: JSON lines on stderr (and the `analyze` report as JSON on stdout)
    #[arg(long, global = true)]
    pub json: bool,

    /// Only report warnings, errors and the final result
    #[arg(short, long, global = true)]
    pub quiet: bool,

    /// Number of worker threads
    #[arg(short = 'w', long, global = true, default_value_t = 3, value_parser = clap::value_parser!(u16).range(1..))]
    pub threads: u16,

    #[command(subcommand)]
    pub command: Command,
}

#[derive(Subcommand, Debug)]
pub enum Command {
    /// Trim Illumina-style short reads (single-end, or paired-end with --r2)
    Short(Box<ShortArgs>),
    /// Trim long reads (Oxford Nanopore, PacBio)
    Long(Box<LongArgs>),
    /// Inspect input files and suggest settings for `simpletrim short` or `simpletrim long`
    Analyze(AnalyzeArgs),
}

/// What to do with a pair when only one mate passes.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum OrphanMode {
    /// Drop the whole pair
    Discard,
    /// Keep the pair; the failing mate is written untrimmed
    KeepBoth,
    /// Write the passing mate to --orphan1/--orphan2
    Orphan,
}

impl OrphanMode {
    pub fn name(self) -> &'static str {
        match self {
            OrphanMode::Discard => "discard",
            OrphanMode::KeepBoth => "keep-both",
            OrphanMode::Orphan => "orphan",
        }
    }
}

/// What `long` does with low-quality stretches inside reads.
#[derive(Clone, Copy, Debug, PartialEq, Eq, ValueEnum)]
pub enum LowQualityMode {
    /// Remove the stretches; every remaining piece becomes a read of its own
    Split,
    /// Remove the stretches and keep only the longest remaining piece
    Longest,
}

impl LowQualityMode {
    pub fn name(self) -> &'static str {
        match self {
            LowQualityMode::Split => "split",
            LowQualityMode::Longest => "longest",
        }
    }
}

fn parse_percent(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("'{s}' is not a number"))?;
    if (0.0..=100.0).contains(&v) {
        Ok(v)
    } else {
        Err(format!("{v} is not between 0 and 100"))
    }
}

fn parse_rate(s: &str) -> Result<f64, String> {
    let v: f64 = s.parse().map_err(|_| format!("'{s}' is not a number"))?;
    if (0.0..=1.0).contains(&v) {
        Ok(v)
    } else {
        Err(format!("{v} is not between 0 and 1"))
    }
}

fn parse_preset(s: &str) -> Result<String, String> {
    if ADAPTER_PRESETS
        .iter()
        .any(|(n, ..)| n.eq_ignore_ascii_case(s))
    {
        Ok(s.to_ascii_lowercase())
    } else {
        let names: Vec<&str> = ADAPTER_PRESETS.iter().map(|(n, ..)| *n).collect();
        Err(format!("unknown preset; available: {}", names.join(", ")))
    }
}

#[derive(Args, Debug)]
pub struct ShortArgs {
    /// Read 1 input (FASTQ, gzip-compressed if the name ends in .gz)
    #[arg(long, value_name = "FILE", help_heading = "Input/output")]
    pub r1: PathBuf,

    /// Read 2 input; enables paired-end mode
    #[arg(
        long,
        value_name = "FILE",
        help_heading = "Input/output",
        requires = "o2"
    )]
    pub r2: Option<PathBuf>,

    /// Read 1 output (gzip-compressed if the name ends in .gz)
    #[arg(long, value_name = "FILE", help_heading = "Input/output")]
    pub o1: PathBuf,

    /// Read 2 output (paired-end only)
    #[arg(
        long,
        value_name = "FILE",
        help_heading = "Input/output",
        requires = "r2"
    )]
    pub o2: Option<PathBuf>,

    /// Paired-end: what happens to a pair when only one read passes
    #[arg(long, value_enum, default_value_t = OrphanMode::Discard, help_heading = "Input/output")]
    pub orphan_mode: OrphanMode,

    /// Output for orphaned reads (--orphan-mode orphan); also gets read 2 orphans unless --orphan2 is set
    #[arg(long, value_name = "FILE", help_heading = "Input/output")]
    pub orphan1: Option<PathBuf>,

    /// Separate output for orphaned read 2 reads (--orphan-mode orphan)
    #[arg(
        long,
        value_name = "FILE",
        help_heading = "Input/output",
        requires = "orphan1"
    )]
    pub orphan2: Option<PathBuf>,

    /// Write the final summary as a JSON file
    #[arg(long, value_name = "FILE", help_heading = "Input/output")]
    pub report: Option<PathBuf>,

    /// Overwrite existing output files
    #[arg(long, help_heading = "Input/output")]
    pub force: bool,

    /// gzip compression level for .gz outputs (1 = fastest, 9 = smallest)
    #[arg(long, default_value_t = 4, value_parser = clap::value_parser!(u32).range(1..=9), help_heading = "Input/output")]
    pub compression_level: u32,

    /// Only process the first N reads (pairs); 0 processes everything
    #[arg(
        long,
        value_name = "N",
        default_value_t = 0,
        help_heading = "Input/output"
    )]
    pub max_reads: u64,

    /// Do not check that read 1 and read 2 names match (by default a mismatch aborts the run)
    #[arg(long, help_heading = "Input/output")]
    pub no_name_check: bool,

    /// Seconds between progress reports; 0 disables them
    #[arg(
        long,
        value_name = "SECONDS",
        default_value_t = 5.0,
        help_heading = "Input/output"
    )]
    pub progress_interval: f64,

    /// Bases to trim from the 5' end of read 1
    #[arg(
        long,
        value_name = "N",
        default_value_t = 0,
        help_heading = "Fixed trimming"
    )]
    pub trim_front1: usize,

    /// Bases to trim from the 3' end of read 1
    #[arg(
        long,
        value_name = "N",
        default_value_t = 0,
        help_heading = "Fixed trimming"
    )]
    pub trim_tail1: usize,

    /// Bases to trim from the 5' end of read 2 [default: --trim-front1]
    #[arg(long, value_name = "N", help_heading = "Fixed trimming")]
    pub trim_front2: Option<usize>,

    /// Bases to trim from the 3' end of read 2 [default: --trim-tail1]
    #[arg(long, value_name = "N", help_heading = "Fixed trimming")]
    pub trim_tail2: Option<usize>,

    /// Slide a window from the 5' end and drop bases while the window mean quality is too low
    #[arg(long, help_heading = "Sliding window quality trimming")]
    pub cut_front: bool,

    /// Slide a window from the 3' end and drop bases while the window mean quality is too low
    #[arg(long, help_heading = "Sliding window quality trimming")]
    pub cut_tail: bool,

    /// Slide a window from the 5' end; at the first low-quality window drop it and everything after (overrides --cut-tail)
    #[arg(long, help_heading = "Sliding window quality trimming")]
    pub cut_right: bool,

    /// Window size for all cut modes
    #[arg(long, value_name = "N", default_value_t = 4, value_parser = clap::value_parser!(u16).range(1..=1000), help_heading = "Sliding window quality trimming")]
    pub cut_window_size: u16,

    /// Required mean phred quality for all cut modes
    #[arg(long, value_name = "Q", default_value_t = 20, value_parser = clap::value_parser!(u8).range(1..=36), help_heading = "Sliding window quality trimming")]
    pub cut_mean_quality: u8,

    /// Window size for --cut-front [default: --cut-window-size]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u16).range(1..=1000), help_heading = "Sliding window quality trimming")]
    pub cut_front_window_size: Option<u16>,

    /// Mean quality for --cut-front [default: --cut-mean-quality]
    #[arg(long, value_name = "Q", value_parser = clap::value_parser!(u8).range(1..=36), help_heading = "Sliding window quality trimming")]
    pub cut_front_mean_quality: Option<u8>,

    /// Window size for --cut-tail [default: --cut-window-size]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u16).range(1..=1000), help_heading = "Sliding window quality trimming")]
    pub cut_tail_window_size: Option<u16>,

    /// Mean quality for --cut-tail [default: --cut-mean-quality]
    #[arg(long, value_name = "Q", value_parser = clap::value_parser!(u8).range(1..=36), help_heading = "Sliding window quality trimming")]
    pub cut_tail_mean_quality: Option<u8>,

    /// Window size for --cut-right [default: --cut-window-size]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u16).range(1..=1000), help_heading = "Sliding window quality trimming")]
    pub cut_right_window_size: Option<u16>,

    /// Mean quality for --cut-right [default: --cut-mean-quality]
    #[arg(long, value_name = "Q", value_parser = clap::value_parser!(u8).range(1..=36), help_heading = "Sliding window quality trimming")]
    pub cut_right_mean_quality: Option<u8>,

    /// Always trim poly-G tails (by default only for two-colour instruments: NextSeq, NovaSeq, ...)
    #[arg(
        long,
        conflicts_with = "disable_trim_poly_g",
        help_heading = "Poly-G trimming"
    )]
    pub trim_poly_g: bool,

    /// Never trim poly-G tails
    #[arg(long, help_heading = "Poly-G trimming")]
    pub disable_trim_poly_g: bool,

    /// Minimum poly-G tail length
    #[arg(long, value_name = "N", default_value_t = 10, value_parser = clap::value_parser!(u16).range(1..), help_heading = "Poly-G trimming")]
    pub poly_g_min_len: u16,

    /// Trim homopolymer tails of any base from the 3' end (after adapter trimming)
    #[arg(long, help_heading = "Poly-X trimming")]
    pub trim_poly_x: bool,

    /// Minimum poly-X tail length
    #[arg(long, value_name = "N", default_value_t = 10, value_parser = clap::value_parser!(u16).range(1..), help_heading = "Poly-X trimming")]
    pub poly_x_min_len: u16,

    /// Disable adapter trimming (and the adapter dimer check)
    #[arg(long, help_heading = "Adapter trimming")]
    pub disable_adapter_trimming: bool,

    /// Read 1 adapter: a sequence, `auto` (detect; for paired-end only with --detect-adapter-pe) or `none`
    #[arg(
        long,
        value_name = "SEQ|auto|none",
        default_value = "auto",
        help_heading = "Adapter trimming"
    )]
    pub adapter_r1: String,

    /// Read 2 adapter: a sequence, `auto` or `none` [default: the read 1 adapter, unless detecting]
    #[arg(
        long,
        value_name = "SEQ|auto|none",
        default_value = "auto",
        help_heading = "Adapter trimming"
    )]
    pub adapter_r2: String,

    /// Use a common adapter pair: truseq, nextera, smallrna, mgi (explicit --adapter-r1/--adapter-r2 take precedence)
    #[arg(long, value_name = "NAME", value_parser = parse_preset, help_heading = "Adapter trimming")]
    pub adapter_preset: Option<String>,

    /// FASTA file with additional adapter sequences, all trimmed from both reads
    #[arg(long, value_name = "FILE", help_heading = "Adapter trimming")]
    pub adapter_fasta: Option<PathBuf>,

    /// Paired-end: also auto-detect adapter sequences, used when reads do not overlap
    #[arg(long, help_heading = "Adapter trimming")]
    pub detect_adapter_pe: bool,

    /// Minimum overlap of an adapter with the read end [default: 4; 5/6 for FASTA lists > 16/256]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u16).range(1..), help_heading = "Adapter trimming")]
    pub adapter_min_overlap: Option<u16>,

    /// Allowed mismatches per compared adapter base
    #[arg(long, value_name = "RATE", default_value_t = 0.125, value_parser = parse_rate, help_heading = "Adapter trimming")]
    pub adapter_error_rate: f64,

    /// Maximum mismatches in an adapter match, regardless of length [default: no limit]
    #[arg(long, value_name = "N", help_heading = "Adapter trimming")]
    pub adapter_max_mismatches: Option<u16>,

    /// Paired-end: minimum overlap of read 1 and read 2
    #[arg(
        long,
        value_name = "N",
        default_value_t = 30,
        help_heading = "Adapter trimming"
    )]
    pub overlap_min_len: u16,

    /// Paired-end: maximum mismatches in the read overlap
    #[arg(
        long,
        value_name = "N",
        default_value_t = 5,
        help_heading = "Adapter trimming"
    )]
    pub overlap_max_mismatches: u16,

    /// Paired-end: maximum mismatch fraction of the read overlap
    #[arg(long, value_name = "RATE", default_value_t = 0.2, value_parser = parse_rate, help_heading = "Adapter trimming")]
    pub overlap_error_rate: f64,

    /// Paired-end: allow one insertion/deletion in the read overlap
    #[arg(long, help_heading = "Adapter trimming")]
    pub allow_gap_overlap_trimming: bool,

    /// Disable the adapter dimer check
    #[arg(long, help_heading = "Adapter dimer check")]
    pub disable_dimer_check: bool,

    /// Reads at most this long after adapter trimming are adapter dimers (all reads of a pair must qualify)
    #[arg(
        long,
        value_name = "N",
        default_value_t = 2,
        help_heading = "Adapter dimer check"
    )]
    pub dimer_max_len: usize,

    /// Disable the read quality filter
    #[arg(long, help_heading = "Filters")]
    pub disable_quality_filter: bool,

    /// Phred quality from which a base counts as qualified
    #[arg(long, value_name = "Q", default_value_t = 15, value_parser = clap::value_parser!(u8).range(0..=93), help_heading = "Filters")]
    pub qualified_quality: u8,

    /// Maximum percentage of unqualified bases
    #[arg(long, value_name = "PERCENT", default_value_t = 40, value_parser = clap::value_parser!(u32).range(0..=100), help_heading = "Filters")]
    pub unqualified_percent_limit: u32,

    /// Maximum number of N bases
    #[arg(long, value_name = "N", default_value_t = 5, help_heading = "Filters")]
    pub n_base_limit: usize,

    /// Minimum mean phred quality; 0 disables this check
    #[arg(long, value_name = "Q", default_value_t = 0, help_heading = "Filters")]
    pub average_quality: u32,

    /// Disable the length filter
    #[arg(long, help_heading = "Filters")]
    pub disable_length_filter: bool,

    /// Discard reads shorter than this after trimming
    #[arg(long, value_name = "N", default_value_t = 15, help_heading = "Filters")]
    pub min_length: usize,

    /// Discard reads longer than this after trimming; 0 means no limit
    #[arg(long, value_name = "N", default_value_t = 0, help_heading = "Filters")]
    pub max_length: usize,
}

/// Options of `long`. Defaults follow fastplong; window and read means of quality values
/// average error probabilities.
#[derive(Args, Debug)]
pub struct LongArgs {
    /// Input (FASTQ, gzip-compressed if the name ends in .gz)
    #[arg(short = 'i', long, value_name = "FILE", help_heading = "Input/output")]
    pub input: PathBuf,

    /// Output (gzip-compressed if the name ends in .gz)
    #[arg(short = 'o', long, value_name = "FILE", help_heading = "Input/output")]
    pub output: PathBuf,

    /// Also write reads and read pieces that fail, with the reason appended to the read name
    #[arg(long, value_name = "FILE", help_heading = "Input/output")]
    pub failed_output: Option<PathBuf>,

    /// Write the final summary as a JSON file
    #[arg(long, value_name = "FILE", help_heading = "Input/output")]
    pub report: Option<PathBuf>,

    /// Overwrite existing output files
    #[arg(long, help_heading = "Input/output")]
    pub force: bool,

    /// gzip compression level for .gz outputs (1 = fastest, 9 = smallest)
    #[arg(long, default_value_t = 4, value_parser = clap::value_parser!(u32).range(1..=9), help_heading = "Input/output")]
    pub compression_level: u32,

    /// Only process the first N reads; 0 processes everything
    #[arg(
        long,
        value_name = "N",
        default_value_t = 0,
        help_heading = "Input/output"
    )]
    pub max_reads: u64,

    /// Seconds between progress reports; 0 disables them
    #[arg(
        long,
        value_name = "SECONDS",
        default_value_t = 5.0,
        help_heading = "Input/output"
    )]
    pub progress_interval: f64,

    /// Bases to trim from the 5' end
    #[arg(
        long,
        value_name = "N",
        default_value_t = 0,
        help_heading = "Fixed trimming"
    )]
    pub trim_front: usize,

    /// Bases to trim from the 3' end
    #[arg(
        long,
        value_name = "N",
        default_value_t = 0,
        help_heading = "Fixed trimming"
    )]
    pub trim_tail: usize,

    /// Slide a window from the 5' end and drop bases while the window mean quality is too low
    #[arg(long, help_heading = "Sliding window quality trimming")]
    pub cut_front: bool,

    /// Slide a window from the 3' end and drop bases while the window mean quality is too low
    #[arg(long, help_heading = "Sliding window quality trimming")]
    pub cut_tail: bool,

    /// Window size for both cut modes
    #[arg(long, value_name = "N", default_value_t = 4, value_parser = clap::value_parser!(u16).range(1..=1000), help_heading = "Sliding window quality trimming")]
    pub cut_window_size: u16,

    /// Required mean quality for both cut modes (error-probability mean; Q17 corresponds to
    /// fastplong's Q20 for plain means)
    #[arg(long, value_name = "Q", default_value_t = 17, value_parser = clap::value_parser!(u8).range(1..=36), help_heading = "Sliding window quality trimming")]
    pub cut_mean_quality: u8,

    /// Window size for --cut-front [default: --cut-window-size]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u16).range(1..=1000), help_heading = "Sliding window quality trimming")]
    pub cut_front_window_size: Option<u16>,

    /// Mean quality for --cut-front [default: --cut-mean-quality]
    #[arg(long, value_name = "Q", value_parser = clap::value_parser!(u8).range(1..=36), help_heading = "Sliding window quality trimming")]
    pub cut_front_mean_quality: Option<u8>,

    /// Window size for --cut-tail [default: --cut-window-size]
    #[arg(long, value_name = "N", value_parser = clap::value_parser!(u16).range(1..=1000), help_heading = "Sliding window quality trimming")]
    pub cut_tail_window_size: Option<u16>,

    /// Mean quality for --cut-tail [default: --cut-mean-quality]
    #[arg(long, value_name = "Q", value_parser = clap::value_parser!(u8).range(1..=36), help_heading = "Sliding window quality trimming")]
    pub cut_tail_mean_quality: Option<u8>,

    /// Trim homopolymer tails of any base from the 3' end (before adapter trimming, as fastplong)
    #[arg(long, help_heading = "Poly-X trimming")]
    pub trim_poly_x: bool,

    /// Minimum poly-X tail length
    #[arg(long, value_name = "N", default_value_t = 10, value_parser = clap::value_parser!(u16).range(1..), help_heading = "Poly-X trimming")]
    pub poly_x_min_len: u16,

    /// Disable adapter trimming (also disables chimera splitting)
    #[arg(long, help_heading = "Adapter trimming")]
    pub disable_adapter_trimming: bool,

    /// Adapter at the read start (5'): a sequence, `auto` (detect) or `none`
    #[arg(
        short = 's',
        long,
        value_name = "SEQ|auto|none",
        default_value = "auto",
        help_heading = "Adapter trimming"
    )]
    pub start_adapter: String,

    /// Adapter at the read end (3'): a sequence, `auto` or `none` [default: the reverse
    /// complement of a given --start-adapter, otherwise detected]
    #[arg(
        short = 'e',
        long,
        value_name = "SEQ|auto|none",
        default_value = "auto",
        help_heading = "Adapter trimming"
    )]
    pub end_adapter: String,

    /// FASTA file with additional sequences, trimmed at both read ends and searched inside reads
    #[arg(long, value_name = "FILE", help_heading = "Adapter trimming")]
    pub adapter_fasta: Option<PathBuf>,

    /// Allowed edits (substitutions, insertions, deletions) per adapter base at the read ends
    #[arg(long, value_name = "RATE", default_value_t = 0.25, value_parser = parse_rate, help_heading = "Adapter trimming")]
    pub distance_threshold: f64,

    /// Minimum adapter overlap with the read when the read starts or ends inside the adapter
    #[arg(long, value_name = "N", default_value_t = 16, value_parser = clap::value_parser!(u16).range(1..), help_heading = "Adapter trimming")]
    pub adapter_min_overlap: u16,

    /// Bases searched for adapters at each read end
    #[arg(long, value_name = "N", default_value_t = 200, value_parser = clap::value_parser!(u32).range(10..), help_heading = "Adapter trimming")]
    pub adapter_window: u32,

    /// Extra bases removed next to every adapter found, for a clean cut
    #[arg(long, value_name = "N", default_value_t = 10, value_parser = clap::value_parser!(u16).range(0..=100), help_heading = "Adapter trimming")]
    pub trimming_extension: u16,

    /// Do not split reads at adapters found inside them
    #[arg(long, help_heading = "Chimeric reads")]
    pub disable_chimera_splitting: bool,

    /// Discard reads with an adapter inside instead of splitting them
    #[arg(
        long,
        conflicts_with = "disable_chimera_splitting",
        help_heading = "Chimeric reads"
    )]
    pub discard_chimeric_reads: bool,

    /// Minimum length of a piece kept after splitting at internal adapters
    #[arg(long, value_name = "N", default_value_t = 20, value_parser = clap::value_parser!(u32).range(1..), help_heading = "Chimeric reads")]
    pub chimera_min_segment_length: u32,

    /// Allowed edits per adapter base for adapters inside reads (stricter than at the ends,
    /// because the whole read is searched)
    #[arg(long, value_name = "RATE", default_value_t = 0.15, value_parser = parse_rate, help_heading = "Chimeric reads")]
    pub chimera_distance_threshold: f64,

    /// Remove low-quality stretches inside reads: `split` keeps every remaining piece as a read
    /// of its own (named <id>_1, <id>_2, ...), `longest` keeps only the longest piece [default: off]
    #[arg(
        long,
        value_enum,
        value_name = "MODE",
        help_heading = "Low-quality stretches"
    )]
    pub low_quality_mode: Option<LowQualityMode>,

    /// Window size for finding low-quality stretches
    #[arg(long, value_name = "N", default_value_t = 100, value_parser = clap::value_parser!(u32).range(5..=1_000_000), help_heading = "Low-quality stretches")]
    pub low_quality_window: u32,

    /// Windows with a mean quality below this are low quality (error-probability mean; Q7
    /// corresponds to fastplong's Q10 for plain means)
    #[arg(long, value_name = "Q", default_value_t = 7, value_parser = clap::value_parser!(u8).range(1..=30), help_heading = "Low-quality stretches")]
    pub low_quality_threshold: u8,

    /// Disable the read quality filter
    #[arg(long, help_heading = "Filters")]
    pub disable_quality_filter: bool,

    /// Phred quality from which a base counts as qualified
    #[arg(long, value_name = "Q", default_value_t = 15, value_parser = clap::value_parser!(u8).range(0..=93), help_heading = "Filters")]
    pub qualified_quality: u8,

    /// Maximum percentage of unqualified bases
    #[arg(long, value_name = "PERCENT", default_value_t = 40, value_parser = clap::value_parser!(u32).range(0..=100), help_heading = "Filters")]
    pub unqualified_percent_limit: u32,

    /// Maximum percentage of N bases
    #[arg(long, value_name = "PERCENT", default_value_t = 10, value_parser = clap::value_parser!(u32).range(0..=100), help_heading = "Filters")]
    pub n_percent_limit: u32,

    /// Maximum number of N bases [default: no limit]
    #[arg(long, value_name = "N", help_heading = "Filters")]
    pub n_base_limit: Option<usize>,

    /// Minimum mean read quality (error-probability mean); 0 disables this check
    #[arg(long, value_name = "Q", default_value_t = 0, value_parser = clap::value_parser!(u8).range(0..=93), help_heading = "Filters")]
    pub average_quality: u8,

    /// Minimum GC content in percent [default: no limit]
    #[arg(long, value_name = "PERCENT", value_parser = parse_percent, help_heading = "Filters")]
    pub min_gc: Option<f64>,

    /// Maximum GC content in percent [default: no limit]
    #[arg(long, value_name = "PERCENT", value_parser = parse_percent, help_heading = "Filters")]
    pub max_gc: Option<f64>,

    /// Discard low-complexity reads (too few bases that differ from the next base)
    #[arg(long, help_heading = "Filters")]
    pub low_complexity_filter: bool,

    /// Minimum percentage of bases that differ from the next base
    #[arg(long, value_name = "PERCENT", default_value_t = 30, value_parser = clap::value_parser!(u32).range(0..=100), help_heading = "Filters")]
    pub complexity_threshold: u32,

    /// Disable the length filter
    #[arg(long, help_heading = "Filters")]
    pub disable_length_filter: bool,

    /// Discard reads (pieces) shorter than this after trimming
    #[arg(long, value_name = "N", default_value_t = 20, help_heading = "Filters")]
    pub min_length: usize,

    /// Discard reads (pieces) longer than this after trimming; 0 means no limit
    #[arg(long, value_name = "N", default_value_t = 0, help_heading = "Filters")]
    pub max_length: usize,
}

#[derive(Args, Debug)]
pub struct AnalyzeArgs {
    /// Read 1 (or single-end) input
    #[arg(long, value_name = "FILE")]
    pub r1: PathBuf,

    /// Read 2 input
    #[arg(long, value_name = "FILE")]
    pub r2: Option<PathBuf>,

    /// Number of reads (pairs) to inspect from the start of the input; 0 reads everything
    #[arg(long, value_name = "N", default_value_t = 1_000_000)]
    pub max_reads: u64,

    /// Long reads: stop after this many bases; 0 means no limit
    #[arg(long, value_name = "N", default_value_t = 500_000_000)]
    pub max_bases: u64,
}

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_definition_is_valid() {
        Cli::command().debug_assert();
    }

    #[test]
    fn parses_paired_command() {
        let cli = Cli::try_parse_from([
            "simpletrim",
            "short",
            "--r1",
            "a.fq.gz",
            "--r2",
            "b.fq.gz",
            "--o1",
            "x.fq.gz",
            "--o2",
            "y.fq.gz",
            "--json",
            "--orphan-mode",
            "keep-both",
        ])
        .unwrap();
        assert!(cli.json);
        let Command::Short(args) = cli.command else {
            panic!("expected short");
        };
        assert_eq!(args.orphan_mode, OrphanMode::KeepBoth);
        assert_eq!(args.compression_level, 4);
    }

    #[test]
    fn o2_requires_r2() {
        assert!(
            Cli::try_parse_from(["simpletrim", "short", "--r1", "a", "--o1", "b", "--o2", "c"])
                .is_err()
        );
        assert!(
            Cli::try_parse_from(["simpletrim", "short", "--r1", "a", "--o1", "b", "--r2", "c"])
                .is_err()
        );
    }

    #[test]
    fn parses_long_command() {
        let cli = Cli::try_parse_from([
            "simpletrim",
            "long",
            "-i",
            "in.fq.gz",
            "-o",
            "out.fq.gz",
            "--low-quality-mode",
            "split",
            "-s",
            "ACGTACGTAC",
        ])
        .unwrap();
        let Command::Long(args) = cli.command else {
            panic!("expected long");
        };
        assert_eq!(args.low_quality_mode, Some(LowQualityMode::Split));
        assert_eq!(args.start_adapter, "ACGTACGTAC");
        assert_eq!(args.end_adapter, "auto");
        assert_eq!(args.min_length, 20);
        assert!(
            Cli::try_parse_from([
                "simpletrim",
                "long",
                "-i",
                "a",
                "-o",
                "b",
                "--discard-chimeric-reads",
                "--disable-chimera-splitting"
            ])
            .is_err()
        );
    }
}
