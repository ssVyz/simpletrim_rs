//! Helpers for the long-read integration tests.
#![allow(dead_code)]

pub mod synth;

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

/// A scratch directory that is removed when dropped.
pub struct TempDir(PathBuf);

impl TempDir {
    pub fn new() -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let n = COUNTER.fetch_add(1, Ordering::SeqCst);
        let dir =
            std::env::temp_dir().join(format!("simpletrim-long-test-{}-{n}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }

    pub fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

pub fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_simpletrim"))
        .args(args)
        .output()
        .expect("failed to run simpletrim")
}

pub fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

/// Writes text, gzip-compressed if the name ends in `.gz`.
pub fn write_text(path: &Path, text: &str) {
    if path.extension().is_some_and(|e| e == "gz") {
        let mut enc = GzEncoder::new(fs::File::create(path).unwrap(), Compression::fast());
        enc.write_all(text.as_bytes()).unwrap();
        enc.finish().unwrap();
    } else {
        fs::write(path, text).unwrap();
    }
}

pub fn read_text(path: &Path) -> String {
    let mut text = String::new();
    if path.extension().is_some_and(|e| e == "gz") {
        MultiGzDecoder::new(fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
    } else {
        text = fs::read_to_string(path).unwrap();
    }
    text
}

#[derive(Clone, Debug, PartialEq)]
pub struct Fq {
    /// Header line without `@`.
    pub header: String,
    pub seq: String,
    pub qual: String,
}

impl Fq {
    /// Read id: the header up to the first whitespace.
    pub fn id(&self) -> &str {
        self.header.split_ascii_whitespace().next().unwrap_or("")
    }
}

pub fn parse_fastq(text: &str) -> Vec<Fq> {
    let lines: Vec<&str> = text.lines().collect();
    lines
        .chunks(4)
        .map(|c| Fq {
            header: c[0].trim_start_matches('@').to_string(),
            seq: c[1].to_string(),
            qual: c[3].to_string(),
        })
        .collect()
}

/// JSON events written to stderr with `--json`.
pub fn events(out: &Output) -> Vec<serde_json::Value> {
    String::from_utf8_lossy(&out.stderr)
        .lines()
        .map(|l| serde_json::from_str(l).expect("JSON event"))
        .collect()
}
