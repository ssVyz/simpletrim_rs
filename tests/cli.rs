//! End-to-end tests running the `simpletrim` binary on small synthetic FASTQ files.

use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Output};
use std::sync::atomic::{AtomicUsize, Ordering};

use flate2::Compression;
use flate2::read::MultiGzDecoder;
use flate2::write::GzEncoder;

const ADAPTER1: &str = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCACGCCACAGGATCTCGTATGCCGTCTTCTGCTTG";
const ADAPTER2: &str = "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGTAGATCTCGGTGGTCGCCGTATCATTAAAAAA";
const GOOD: u8 = b'I';

/// A scratch directory that is removed when dropped.
struct TempDir(PathBuf);

impl TempDir {
    fn new() -> Self {
        static COUNTER: AtomicUsize = AtomicUsize::new(0);
        let n = COUNTER.fetch_add(1, Ordering::SeqCst);
        let dir = std::env::temp_dir().join(format!("simpletrim-test-{}-{n}", std::process::id()));
        let _ = fs::remove_dir_all(&dir);
        fs::create_dir_all(&dir).unwrap();
        TempDir(dir)
    }

    fn path(&self, name: &str) -> PathBuf {
        self.0.join(name)
    }
}

impl Drop for TempDir {
    fn drop(&mut self) {
        let _ = fs::remove_dir_all(&self.0);
    }
}

#[derive(Clone, Debug, PartialEq)]
struct Rec {
    name: String,
    seq: String,
    qual: String,
}

fn rec(name: &str, seq: &str, qual: Option<String>) -> Rec {
    Rec {
        name: name.to_string(),
        seq: seq.to_string(),
        qual: qual.unwrap_or_else(|| String::from_utf8(vec![GOOD; seq.len()]).unwrap()),
    }
}

fn random_seq(state: &mut u64, len: usize) -> String {
    (0..len)
        .map(|_| {
            *state = state
                .wrapping_mul(6364136223846793005)
                .wrapping_add(1442695040888963407);
            b"ACGT"[(*state >> 33) as usize % 4] as char
        })
        .collect()
}

fn revcomp(s: &str) -> String {
    s.chars()
        .rev()
        .map(|c| match c {
            'A' => 'T',
            'T' => 'A',
            'C' => 'G',
            'G' => 'C',
            _ => 'N',
        })
        .collect()
}

fn fastq_text(recs: &[Rec]) -> String {
    recs.iter()
        .map(|r| format!("@{}\n{}\n+\n{}\n", r.name, r.seq, r.qual))
        .collect()
}

fn write_fastq(path: &Path, recs: &[Rec]) {
    let text = fastq_text(recs);
    if path.extension().is_some_and(|e| e == "gz") {
        let mut enc = GzEncoder::new(fs::File::create(path).unwrap(), Compression::fast());
        enc.write_all(text.as_bytes()).unwrap();
        enc.finish().unwrap();
    } else {
        fs::write(path, text).unwrap();
    }
}

fn read_fastq(path: &Path) -> Vec<Rec> {
    let mut text = String::new();
    if path.extension().is_some_and(|e| e == "gz") {
        MultiGzDecoder::new(fs::File::open(path).unwrap())
            .read_to_string(&mut text)
            .unwrap();
    } else {
        text = fs::read_to_string(path).unwrap();
    }
    let lines: Vec<&str> = text.lines().collect();
    lines
        .chunks(4)
        .map(|c| Rec {
            name: c[0].trim_start_matches('@').to_string(),
            seq: c[1].to_string(),
            qual: c[3].to_string(),
        })
        .collect()
}

fn run(args: &[&str]) -> Output {
    Command::new(env!("CARGO_BIN_EXE_simpletrim"))
        .args(args)
        .output()
        .expect("failed to run simpletrim")
}

fn s(p: &Path) -> &str {
    p.to_str().unwrap()
}

/// Five read pairs covering the main cases.
fn paired_input(dir: &TempDir, gz: bool) -> (PathBuf, PathBuf) {
    let mut st = 99;
    let ext = if gz { "fastq.gz" } else { "fastq" };
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    // p1: long insert, no adapter
    let ins = random_seq(&mut st, 200);
    r1.push(rec("p1/1", &ins[..150], None));
    r2.push(rec("p1/2", &revcomp(&ins)[..150], None));
    // p2: 60 bp insert, both reads run into adapter
    let ins = random_seq(&mut st, 60);
    r1.push(rec("p2/1", &format!("{ins}{}", &ADAPTER1[..60]), None));
    r2.push(rec(
        "p2/2",
        &format!("{}{}", revcomp(&ins), &ADAPTER2[..60]),
        None,
    ));
    // p3: read 2 has low quality everywhere
    let ins = random_seq(&mut st, 200);
    r1.push(rec("p3/1", &ins[..150], None));
    r2.push(rec("p3/2", &revcomp(&ins)[..150], Some("#".repeat(150))));
    // p4: read 1 has too many N
    let ins = random_seq(&mut st, 200);
    let mut s1 = ins[..150].to_string();
    s1.replace_range(10..20, "NNNNNNNNNN");
    r1.push(rec("p4/1", &s1, None));
    r2.push(rec("p4/2", &revcomp(&ins)[..150], None));
    // p5: adapter dimer (no insert)
    r1.push(rec(
        "p5/1",
        &format!("{ADAPTER1}{}", random_seq(&mut st, 90))[..150],
        None,
    ));
    r2.push(rec(
        "p5/2",
        &format!("{ADAPTER2}{}", random_seq(&mut st, 90))[..150],
        None,
    ));
    let (p1, p2) = (
        dir.path(&format!("in_R1.{ext}")),
        dir.path(&format!("in_R2.{ext}")),
    );
    write_fastq(&p1, &r1);
    write_fastq(&p2, &r2);
    (p1, p2)
}

fn names(recs: &[Rec]) -> Vec<&str> {
    recs.iter().map(|r| r.name.as_str()).collect()
}

const ADAPTER_ARGS: [&str; 4] = [
    "--adapter-r1",
    "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA",
    "--adapter-r2",
    "AGATCGGAAGAGCGTCGTGTAGGGAAAGAGTGT",
];

#[test]
fn paired_discard_mode() {
    let dir = TempDir::new();
    let (r1, r2) = paired_input(&dir, true);
    let (o1, o2) = (dir.path("o1.fastq.gz"), dir.path("o2.fastq.gz"));
    let mut args = vec![
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
    ];
    args.extend(ADAPTER_ARGS);
    let out = run(&args);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let (a, b) = (read_fastq(&o1), read_fastq(&o2));
    assert_eq!(names(&a), ["p1/1", "p2/1"]);
    assert_eq!(names(&b), ["p1/2", "p2/2"]);
    assert_eq!(a[0].seq.len(), 150);
    // the 60 bp insert is recovered from both reads by overlap analysis
    assert_eq!(a[1].seq.len(), 60);
    assert_eq!(b[1].seq.len(), 60);
    assert_eq!(a[1].seq, revcomp(&b[1].seq));
}

#[test]
fn paired_keep_both_mode() {
    let dir = TempDir::new();
    let (r1, r2) = paired_input(&dir, false);
    let (o1, o2) = (dir.path("o1.fastq"), dir.path("o2.fastq"));
    let mut args = vec![
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
        "--orphan-mode",
        "keep-both",
    ];
    args.extend(ADAPTER_ARGS);
    let out = run(&args);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let (a, b) = (read_fastq(&o1), read_fastq(&o2));
    // the dimer pair (p5) is dropped even in keep-both mode
    assert_eq!(names(&a), ["p1/1", "p2/1", "p3/1", "p4/1"]);
    assert_eq!(names(&b), ["p1/2", "p2/2", "p3/2", "p4/2"]);
    // failing mates are written untrimmed
    let input2 = read_fastq(&r2);
    assert_eq!(b[2], input2[2]);
    let input1 = read_fastq(&r1);
    assert_eq!(a[3], input1[3]);
}

#[test]
fn paired_orphan_mode() {
    let dir = TempDir::new();
    let (r1, r2) = paired_input(&dir, false);
    let (o1, o2) = (dir.path("o1.fastq"), dir.path("o2.fastq"));
    let (u1, u2) = (dir.path("orphan1.fastq"), dir.path("orphan2.fastq.gz"));
    let mut args = vec![
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
        "--orphan-mode",
        "orphan",
        "--orphan1",
        s(&u1),
        "--orphan2",
        s(&u2),
    ];
    args.extend(ADAPTER_ARGS);
    let out = run(&args);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(names(&read_fastq(&o1)), ["p1/1", "p2/1"]);
    assert_eq!(names(&read_fastq(&u1)), ["p3/1"]);
    assert_eq!(names(&read_fastq(&u2)), ["p4/2"]);
}

#[test]
fn orphans_share_one_file_without_orphan2() {
    let dir = TempDir::new();
    let (r1, r2) = paired_input(&dir, false);
    let (o1, o2, u1) = (
        dir.path("o1.fastq"),
        dir.path("o2.fastq"),
        dir.path("orphans.fastq"),
    );
    let out = run(&[
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
        "--orphan-mode",
        "orphan",
        "--orphan1",
        s(&u1),
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert_eq!(names(&read_fastq(&u1)), ["p3/1", "p4/2"]);
}

#[test]
fn single_end_with_adapter_and_fixed_trimming() {
    let dir = TempDir::new();
    let mut st = 5;
    let ins = random_seq(&mut st, 80);
    let input = dir.path("in.fastq");
    write_fastq(
        &input,
        &[
            rec("a", &format!("{ins}{}", &ADAPTER1[..60]), None),
            rec("b", &random_seq(&mut st, 140), None),
        ],
    );
    let o1 = dir.path("out.fastq");
    let out = run(&[
        "short",
        "--r1",
        s(&input),
        "--o1",
        s(&o1),
        "--adapter-preset",
        "truseq",
        "--trim-front1",
        "5",
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let recs = read_fastq(&o1);
    assert_eq!(recs[0].seq, ins[5..]);
    assert_eq!(recs[1].seq.len(), 135);
}

#[test]
fn single_end_auto_detects_adapter() {
    let dir = TempDir::new();
    let mut st = 11;
    let recs: Vec<Rec> = (0..12_000)
        .map(|i| {
            let len = if i % 3 == 0 { 50 + i % 60 } else { 150 };
            let seq = format!(
                "{}{ADAPTER1}{}",
                random_seq(&mut st, len),
                random_seq(&mut st, 100)
            );
            rec(&format!("r{i}"), &seq[..150], None)
        })
        .collect();
    let input = dir.path("in.fastq.gz");
    write_fastq(&input, &recs);
    let o1 = dir.path("out.fastq.gz");
    let out = run(&["short", "--r1", s(&input), "--o1", s(&o1), "--json"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    let adapter_event = stderr
        .lines()
        .map(|l| serde_json::from_str::<serde_json::Value>(l).unwrap())
        .find(|v| v["event"] == "adapter")
        .expect("adapter event");
    assert_eq!(
        adapter_event["sequence"],
        "AGATCGGAAGAGCACACGTCTGAACTCCAGTCA"
    );
    let trimmed = read_fastq(&o1);
    assert_eq!(trimmed[0].seq.len(), 50);
}

#[test]
fn refuses_to_overwrite_without_force() {
    let dir = TempDir::new();
    let (r1, r2) = paired_input(&dir, false);
    let (o1, o2) = (dir.path("o1.fastq"), dir.path("o2.fastq"));
    fs::write(&o1, "keep me").unwrap();
    let out = run(&[
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
    ]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("--force"));
    assert_eq!(fs::read_to_string(&o1).unwrap(), "keep me");
    assert!(!o2.exists());
    let out = run(&[
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
        "--force",
    ]);
    assert!(out.status.success());
    assert_ne!(fs::read_to_string(&o1).unwrap(), "keep me");
}

#[test]
fn name_mismatch_aborts_and_removes_outputs() {
    let dir = TempDir::new();
    let (r1, r2) = (dir.path("a.fastq"), dir.path("b.fastq"));
    write_fastq(&r1, &[rec("x/1", "ACGTACGTACGTACGTACGT", None)]);
    write_fastq(&r2, &[rec("y/2", "ACGTACGTACGTACGTACGT", None)]);
    let (o1, o2) = (dir.path("o1.fastq"), dir.path("o2.fastq"));
    let out = run(&[
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
    ]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("do not match"));
    assert!(!o1.exists() && !o2.exists());
    let out = run(&[
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
        "--no-name-check",
    ]);
    assert!(out.status.success());
}

#[test]
fn different_read_counts_fail() {
    let dir = TempDir::new();
    let (r1, r2) = (dir.path("a.fastq"), dir.path("b.fastq"));
    write_fastq(
        &r1,
        &[
            rec("x/1", "ACGTACGTACGTACGTACGT", None),
            rec("z/1", "ACGTACGTACGTACGTACGT", None),
        ],
    );
    write_fastq(&r2, &[rec("x/2", "ACGTACGTACGTACGTACGT", None)]);
    let (o1, o2) = (dir.path("o1.fastq"), dir.path("o2.fastq"));
    let out = run(&[
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
    ]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("different numbers of reads"));
}

#[test]
fn json_events_and_report() {
    let dir = TempDir::new();
    let (r1, r2) = paired_input(&dir, true);
    let (o1, o2, report) = (
        dir.path("o1.fq.gz"),
        dir.path("o2.fq.gz"),
        dir.path("report.json"),
    );
    let out = run(&[
        "short",
        "--r1",
        s(&r1),
        "--r2",
        s(&r2),
        "--o1",
        s(&o1),
        "--o2",
        s(&o2),
        "--json",
        "--report",
        s(&report),
    ]);
    assert!(out.status.success());
    assert!(out.stdout.is_empty());
    let events: Vec<serde_json::Value> = String::from_utf8_lossy(&out.stderr)
        .lines()
        .map(|l| serde_json::from_str(l).expect("every stderr line is JSON"))
        .collect();
    assert_eq!(events.first().unwrap()["event"], "start");
    let summary = events.last().unwrap();
    assert_eq!(summary["event"], "summary");
    assert_eq!(summary["records_in"], 5);
    assert_eq!(summary["pairs"]["pairs"], 5);
    let file: serde_json::Value =
        serde_json::from_str(&fs::read_to_string(&report).unwrap()).unwrap();
    assert_eq!(file["records_written"], summary["records_written"]);
}

#[test]
fn empty_input_gives_valid_empty_gzip() {
    let dir = TempDir::new();
    let input = dir.path("empty.fastq");
    fs::write(&input, "").unwrap();
    let o1 = dir.path("out.fastq.gz");
    let out = run(&["short", "--r1", s(&input), "--o1", s(&o1)]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    assert!(read_fastq(&o1).is_empty());
}

#[test]
fn usage_errors_exit_with_2() {
    let out = run(&["short", "--r1", "x.fastq"]);
    assert_eq!(out.status.code(), Some(2));
    let out = run(&["short", "--r1", "x.fastq", "--json"]);
    assert_eq!(out.status.code(), Some(2));
    let v: serde_json::Value =
        serde_json::from_str(String::from_utf8_lossy(&out.stderr).trim()).unwrap();
    assert_eq!(v["event"], "error");
}

#[test]
fn long_is_not_implemented() {
    let out = run(&["long", "--input", "x"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("not implemented"));
}

#[test]
fn analyze_reports_json() {
    let dir = TempDir::new();
    let (r1, r2) = paired_input(&dir, true);
    let out = run(&["analyze", "--r1", s(&r1), "--r2", s(&r2), "--json"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(report["paired"], true);
    assert_eq!(report["sample"]["records"], 5);
    assert_eq!(report["pairing"]["name_mismatches"], 0);
    let cmd = report["suggested_command"].as_str().unwrap();
    assert!(cmd.starts_with("simpletrim short --r1"));
    assert!(cmd.contains("in_R1.trimmed.fastq.gz"));
}

#[test]
fn analyze_recognises_nanopore() {
    let dir = TempDir::new();
    let input = dir.path("ont.fastq");
    let mut st = 3;
    write_fastq(
        &input,
        &[rec(
            "0a1b2c3d-1234-4abc-9def-0123456789ab runid=abc read=1 ch=5 start_time=2024-01-01T00:00:00Z",
            &random_seq(&mut st, 2000),
            None,
        )],
    );
    let out = run(&["analyze", "--r1", s(&input), "--json"]);
    assert!(out.status.success());
    let report: serde_json::Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(report["status"], "not_implemented");
    assert_eq!(report["platform"]["platform"], "oxford_nanopore");
}
