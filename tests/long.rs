//! End-to-end tests of `simpletrim long` and of `analyze` on long reads, using synthetic
//! nanopore-like reads with a known structure.

mod common;

use std::collections::HashMap;

use common::synth::{Library, Params, Rng, SimRead, simulate_many, to_fastq};
use common::{Fq, TempDir, events, parse_fastq, read_text, run, s, write_text};
use serde_json::Value;

/// The read an output record comes from (its id without a `_<n>` piece suffix), and where the
/// record lies in that read.
fn locate<'a>(rec: &Fq, reads: &HashMap<&str, &'a SimRead>) -> (&'a SimRead, usize, usize) {
    let id = rec.id();
    let name = match id.rsplit_once('_') {
        Some((n, k)) if !reads.contains_key(id) && k.parse::<usize>().is_ok() => n,
        _ => id,
    };
    let read = reads[name];
    let seq = std::str::from_utf8(&read.seq).unwrap();
    let qual = std::str::from_utf8(&read.qual).unwrap();
    let pos = seq
        .match_indices(rec.seq.as_str())
        .map(|(p, _)| p)
        .find(|&p| qual[p..p + rec.seq.len()] == rec.qual)
        .expect("record is a piece of its read");
    (read, pos, pos + rec.seq.len())
}

fn by_name(reads: &[SimRead]) -> HashMap<&str, &SimRead> {
    reads.iter().map(|r| (r.name.as_str(), r)).collect()
}

#[test]
fn trims_synthetic_reads_with_detected_adapters() {
    let dir = TempDir::new();
    let mut rng = Rng::new(21);
    let lib = Library::random(&mut rng);
    let p = Params {
        chimera_rate: 0.05,
        ..Params::default()
    };
    let reads = simulate_many(22, &lib, &p, 3000);
    let input = dir.path("in.fastq.gz");
    write_text(&input, &to_fastq(&reads));
    let output = dir.path("out.fastq.gz");
    let report = dir.path("report.json");
    let out = run(&[
        "long",
        "-i",
        s(&input),
        "-o",
        s(&output),
        "--report",
        s(&report),
        "--json",
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let ev = events(&out);
    let adapters: Vec<&Value> = ev.iter().filter(|e| e["event"] == "adapter").collect();
    assert_eq!(adapters.len(), 2);
    assert!(adapters.iter().all(|e| e["sequence"].is_string()));
    let summary = ev.last().unwrap();
    assert_eq!(summary["event"], "summary");
    assert_eq!(summary["reads_in"], 3000);
    let rep: Value = serde_json::from_str(&read_text(&report)).unwrap();
    assert_eq!(rep["records_written"], summary["records_written"]);
    assert_eq!(rep["settings"]["chimeras"]["mode"], "split");

    let names = by_name(&reads);
    let recs = parse_fastq(&read_text(&output));
    assert_eq!(
        recs.len() as u64,
        summary["records_written"].as_u64().unwrap()
    );
    let (mut good, mut single) = (0, 0);
    let mut chimera_pieces: HashMap<&str, usize> = HashMap::new();
    for rec in &recs {
        let (read, a, b) = locate(rec, &names);
        if read.inserts.len() == 1 {
            single += 1;
            let (s0, e0) = read.inserts[0];
            // no adapter bases left (3 bp slack), at most 30 insert bases lost per end
            if a + 3 >= s0 && a <= s0 + 30 && b <= e0 + 3 && b + 30 >= e0 {
                good += 1;
            }
        } else {
            *chimera_pieces.entry(read.name.as_str()).or_default() += 1;
        }
    }
    assert!(good as f64 >= 0.97 * single as f64, "{good} of {single}");
    let chimeras = reads.iter().filter(|r| r.inserts.len() == 2).count();
    let split = chimera_pieces.values().filter(|&&n| n == 2).count();
    assert!(
        chimeras > 50 && split as f64 >= 0.8 * chimeras as f64,
        "{split} of {chimeras}"
    );
    assert!(recs.iter().any(|r| r.id().ends_with("_2")));
}

#[test]
fn splits_reads_at_low_quality_stretches() {
    let dir = TempDir::new();
    let mut rng = Rng::new(31);
    let lib = Library::random(&mut rng);
    let p = Params {
        insert_len: (2000, 4000),
        low_quality_rate: 1.0,
        start_adapter_rate: 0.0,
        end_adapter_rate: 0.0,
        ..Params::default()
    };
    let reads = simulate_many(32, &lib, &p, 200);
    let input = dir.path("in.fastq");
    write_text(&input, &to_fastq(&reads));
    let names = by_name(&reads);

    let output = dir.path("split.fastq");
    let out = run(&[
        "long",
        "-i",
        s(&input),
        "-o",
        s(&output),
        "--disable-adapter-trimming",
        "--low-quality-mode",
        "split",
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let recs = parse_fastq(&read_text(&output));
    let mut pieces: HashMap<&str, Vec<(usize, usize)>> = HashMap::new();
    for rec in &recs {
        let (read, a, b) = locate(rec, &names);
        let (ls, le) = read.low_quality[0];
        // pieces stay out of the poor stretch
        assert!(
            b <= ls + 20 || a + 20 >= le,
            "{} {a}-{b} vs {ls}-{le}",
            rec.id()
        );
        pieces.entry(read.name.as_str()).or_default().push((a, b));
    }
    let two = pieces.values().filter(|p| p.len() == 2).count();
    assert!(two >= 190, "{two} of 200 reads split in two");
    let kept: usize = pieces.values().flatten().map(|(a, b)| b - a).sum();
    let total: usize = reads.iter().map(|r| r.seq.len() - 300).sum();
    assert!(kept as f64 > 0.95 * total as f64, "{kept} of {total}");

    let output = dir.path("longest.fastq");
    let out = run(&[
        "long",
        "-i",
        s(&input),
        "-o",
        s(&output),
        "--disable-adapter-trimming",
        "--low-quality-mode",
        "longest",
    ]);
    assert!(out.status.success());
    let recs = parse_fastq(&read_text(&output));
    assert_eq!(recs.len(), 200);
    assert!(recs.iter().all(|r| names.contains_key(r.id())));
}

#[test]
fn failed_output_has_reasons() {
    let dir = TempDir::new();
    let mut rng = Rng::new(41);
    let good = String::from_utf8(rng.seq(500)).unwrap();
    let ns = format!("{}{}", "N".repeat(100), &good[..300]);
    let text = format!(
        "@good x=1\n{good}\n+\n{}\n@short\nACGTACGTAC\n+\nIIIIIIIIII\n@nread\n{ns}\n+\n{}\n",
        "5".repeat(500),
        "5".repeat(400)
    );
    let input = dir.path("in.fastq");
    write_text(&input, &text);
    let (output, failed) = (dir.path("out.fastq"), dir.path("failed.fastq.gz"));
    let out = run(&[
        "long",
        "-i",
        s(&input),
        "-o",
        s(&output),
        "--failed-output",
        s(&failed),
        "--disable-adapter-trimming",
    ]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let ok = parse_fastq(&read_text(&output));
    assert_eq!(ok.len(), 1);
    assert_eq!(ok[0].header, "good x=1");
    let bad: Vec<String> = parse_fastq(&read_text(&failed))
        .into_iter()
        .map(|r| r.header)
        .collect();
    assert_eq!(
        bad,
        ["short failed_too_short", "nread failed_too_many_n_bases"]
    );
}

#[test]
fn refuses_to_overwrite_and_rejects_mixed_rna() {
    let dir = TempDir::new();
    let input = dir.path("in.fastq");
    write_text(
        &input,
        "@r1\nACGTACGTACGTACGTACGTACGT\n+\nIIIIIIIIIIIIIIIIIIIIIIII\n",
    );
    let output = dir.path("out.fastq");
    write_text(&output, "keep me");
    let out = run(&["long", "-i", s(&input), "-o", s(&output)]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("already exists"));
    assert_eq!(read_text(&output), "keep me");

    let mixed = dir.path("mixed.fastq");
    write_text(&mixed, "@r1\nACGTACGU\n+\nIIIIIIII\n");
    let out = run(&["long", "-i", s(&mixed), "-o", s(&dir.path("o2.fastq"))]);
    assert_eq!(out.status.code(), Some(1));
    assert!(String::from_utf8_lossy(&out.stderr).contains("both U and T"));
}

#[test]
fn analyze_reports_long_reads() {
    let dir = TempDir::new();
    let mut rng = Rng::new(51);
    let lib = Library::random(&mut rng);
    let mut reads = simulate_many(52, &lib, &Params::default(), 2000);
    for (i, r) in reads.iter_mut().enumerate() {
        r.name = format!("{i:08x}-0000-4000-8000-{i:012x} runid=run7 flow_cell_id=FAB123 ch=1");
    }
    let input = dir.path("ont.fastq.gz");
    write_text(&input, &to_fastq(&reads));
    let out = run(&["analyze", "--r1", s(&input), "--json"]);
    assert!(
        out.status.success(),
        "{}",
        String::from_utf8_lossy(&out.stderr)
    );
    let report: Value = serde_json::from_slice(&out.stdout).unwrap();
    assert_eq!(report["mode"], "long");
    assert_eq!(report["platform"]["run"]["flow_cell"], "FAB123");
    assert_eq!(report["sample"]["records"], 2000);
    assert!(report["reads"]["length"]["n50"].as_u64().unwrap() > 1000);
    let ad = &report["adapters"];
    assert!(ad["start"]["detected"].is_string());
    assert!(ad["end"]["detected"].is_string());
    assert!(ad["start"]["reads_with_adapter_fraction"].as_f64().unwrap() > 0.8);
    let cmd = report["suggested_command"].as_str().unwrap();
    assert!(cmd.starts_with("simpletrim long -i"), "{cmd}");
    assert!(cmd.contains("--start-adapter"));

    // the text report goes to stdout as well
    let out = run(&["analyze", "--r1", s(&input)]);
    assert!(String::from_utf8_lossy(&out.stdout).contains("Suggested command"));
}

/// Writes a synthetic reference dataset to `examples/` (git-ignored): about 200 Mbp of
/// nanopore-like reads with adapters, chimeras and poor stretches, plus a truth table.
/// Run with `cargo test --release --test long reference_dataset -- --ignored`.
#[test]
#[ignore]
fn reference_dataset() {
    let dir = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("examples");
    std::fs::create_dir_all(&dir).unwrap();
    let mut rng = Rng::new(2026);
    let lib = Library::random(&mut rng);
    let p = Params {
        insert_len: (500, 40_000),
        chimera_rate: 0.03,
        low_quality_rate: 0.1,
        ..Params::default()
    };
    let reads = simulate_many(2027, &lib, &p, 20_000);
    write_text(&dir.join("synthetic_long.fastq.gz"), &to_fastq(&reads));
    let mut truth = format!(
        "# start_adapter={} end_adapter={}\n# coordinates are 0-based, end exclusive\nname\tlength\tinserts\tlow_quality\tstart_adapter\tend_adapter\n",
        String::from_utf8_lossy(&lib.start_adapter),
        String::from_utf8_lossy(&lib.end_adapter)
    );
    let ranges = |v: &[(usize, usize)]| {
        v.iter()
            .map(|(a, b)| format!("{a}-{b}"))
            .collect::<Vec<_>>()
            .join(",")
    };
    for r in &reads {
        truth.push_str(&format!(
            "{}\t{}\t{}\t{}\t{}\t{}\n",
            r.name,
            r.seq.len(),
            ranges(&r.inserts),
            ranges(&r.low_quality),
            r.has_start_adapter,
            r.has_end_adapter
        ));
    }
    std::fs::write(dir.join("synthetic_long.truth.tsv"), truth).unwrap();
}
