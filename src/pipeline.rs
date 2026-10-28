//! The reader → workers → writer pipeline shared by the trimming subcommands.
//!
//! One reader thread per input parses chunks; chunks are processed and (for .gz outputs)
//! compressed in parallel on the worker pool; a writer thread writes the results in input
//! order. The number of chunks in flight is bounded, which bounds memory use.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write;
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};
use std::sync::mpsc::{Receiver, SyncSender, channel, sync_channel};
use std::time::{Duration, Instant};

use serde_json::json;

use crate::Result;
use crate::fastq::{Chunk, feed_chunks};
use crate::io::{Input, OutputFiles, gzip_member, is_gz};
use crate::log::Logger;
use crate::stats::group;

/// Statistics collected per chunk and merged in input order.
pub trait Merge: Default + Send {
    fn merge(&mut self, other: &Self);
}

/// An output file being written.
pub struct Sink {
    file: File,
    path: PathBuf,
    gz: bool,
    written: u64,
}

impl Sink {
    /// Creates the output file; it is registered in `files`, so it is removed if the run fails.
    pub fn create(files: &mut OutputFiles, path: &Path, force: bool) -> Result<Sink> {
        Ok(Sink {
            file: files.create(path, force)?,
            path: path.to_path_buf(),
            gz: is_gz(path),
            written: 0,
        })
    }

    fn write(&mut self, data: &[u8]) -> Result<()> {
        self.file
            .write_all(data)
            .map_err(|e| format!("cannot write {}: {e}", self.path.display()))?;
        self.written += data.len() as u64;
        Ok(())
    }

    fn finish(mut self) -> Result<()> {
        if self.gz && self.written == 0 {
            // an empty file is not valid gzip; write an empty member instead
            self.write(&gzip_member(b"", 1))?;
        }
        self.file
            .flush()
            .map_err(|e| format!("cannot write {}: {e}", self.path.display()))
    }
}

/// How the pipeline reads and reports.
pub struct Options<'a> {
    pub threads: usize,
    /// Stop after this many records (pairs); 0 reads everything.
    pub max_reads: u64,
    pub chunk_records: usize,
    /// Byte limit per chunk. Only used for a single input: paired chunks must hold the same
    /// number of records.
    pub chunk_bytes: usize,
    /// Seconds between progress reports; zero, negative or invalid values disable them.
    pub progress_interval: f64,
    /// Unit in progress messages (`reads`, `pairs`).
    pub unit: &'static str,
    pub compression_level: u32,
    pub log: &'a Logger,
}

struct Done<S> {
    index: usize,
    records: u64,
    outputs: Vec<Vec<u8>>,
    stats: S,
}

/// Periodic progress reports from the writer thread.
struct Progress<'a> {
    log: &'a Logger,
    interval: Option<Duration>,
    start: Instant,
    last: Instant,
    unit: &'static str,
    bytes_read: Arc<AtomicU64>,
    file_size: u64,
}

impl Progress<'_> {
    fn update(&mut self, records: u64) {
        let Some(interval) = self.interval else {
            return;
        };
        if self.last.elapsed() < interval {
            return;
        }
        self.last = Instant::now();
        let secs = self.start.elapsed().as_secs_f64().max(1e-9);
        let rate = records as f64 / secs;
        let fraction = (self.file_size > 0).then(|| {
            (self.bytes_read.load(Ordering::Relaxed) as f64 / self.file_size as f64).min(1.0)
        });
        let pct = fraction
            .map(|f| format!(" ({:.1}% of input)", f * 100.0))
            .unwrap_or_default();
        self.log.progress(
            &format!(
                "processed {} {}{pct}, {} {}/s",
                group(records),
                self.unit,
                group(rate as u64),
                self.unit
            ),
            json!({
                "records": records,
                "unit": self.unit,
                "input_fraction": fraction.map(|f| (f * 10000.0).round() / 10000.0),
                "records_per_s": rate.round(),
            }),
        );
    }
}

fn write_loop<S: Merge>(
    rx: Receiver<Result<Done<S>>>,
    tokens: SyncSender<()>,
    mut sinks: Vec<Option<Sink>>,
    mut progress: Progress,
) -> Result<S> {
    let mut pending: BTreeMap<usize, Done<S>> = BTreeMap::new();
    let mut next = 0;
    let mut stats = S::default();
    let mut records = 0u64;
    for msg in rx {
        let done = msg?;
        pending.insert(done.index, done);
        while let Some(done) = pending.remove(&next) {
            for (slot, data) in done.outputs.iter().enumerate() {
                if !data.is_empty() {
                    sinks[slot]
                        .as_mut()
                        .expect("output slot without file")
                        .write(data)?;
                }
            }
            stats.merge(&done.stats);
            records += done.records;
            next += 1;
            let _ = tokens.send(());
            progress.update(records);
        }
    }
    for sink in sinks.into_iter().flatten() {
        sink.finish()?;
    }
    Ok(stats)
}

fn receive(rx: Option<&Receiver<Result<Chunk>>>) -> Result<Option<Chunk>> {
    match rx.map(|rx| rx.recv()) {
        Some(Ok(chunk)) => chunk.map(Some),
        _ => Ok(None),
    }
}

/// Runs the pipeline over one input, or two inputs read in lockstep (paired-end).
///
/// `process` turns a chunk (and the matching chunk of the second input) into one output
/// buffer per sink, plus the statistics of the chunk. Buffers for `.gz` sinks are compressed
/// on the worker thread.
pub fn run<S, F>(
    opts: &Options,
    in1: Input,
    in2: Option<Input>,
    sinks: Vec<Option<Sink>>,
    process: F,
) -> Result<S>
where
    S: Merge,
    F: Fn(&Chunk, Option<&Chunk>) -> Result<(Vec<Vec<u8>>, S)> + Sync,
{
    let paired = in2.is_some();
    let pool = rayon::ThreadPoolBuilder::new()
        .num_threads(opts.threads)
        .thread_name(|i| format!("simpletrim-worker-{i}"))
        .build()
        .map_err(|e| format!("cannot start worker threads: {e}"))?;
    let max_in_flight = opts.threads * 2 + 2;
    let chunk_bytes = if paired { usize::MAX } else { opts.chunk_bytes };
    let gz: Vec<bool> = sinks
        .iter()
        .map(|s| s.as_ref().is_some_and(|s| s.gz))
        .collect();
    let level = opts.compression_level;
    let progress = Progress {
        log: opts.log,
        interval: Duration::try_from_secs_f64(opts.progress_interval)
            .ok()
            .filter(|d| !d.is_zero()),
        start: Instant::now(),
        last: Instant::now(),
        unit: opts.unit,
        bytes_read: in1.bytes_read.clone(),
        file_size: in1.file_size,
    };
    let (max_reads, chunk_records) = (opts.max_reads, opts.chunk_records);
    let process = &process;
    let gz = &gz;

    std::thread::scope(|scope| {
        let (tx1, rx1) = sync_channel(max_in_flight);
        let r1_thread = scope
            .spawn(move || feed_chunks(in1.reader, tx1, max_reads, chunk_records, chunk_bytes));
        let (rx2, r2_thread) = match in2 {
            Some(in2) => {
                let (tx2, rx2) = sync_channel(max_in_flight);
                let t = scope.spawn(move || {
                    feed_chunks(in2.reader, tx2, max_reads, chunk_records, chunk_bytes)
                });
                (Some(rx2), Some(t))
            }
            None => (None, None),
        };
        let (done_tx, done_rx) = channel::<Result<Done<S>>>();
        let (token_tx, token_rx) = sync_channel::<()>(max_in_flight);
        for _ in 0..max_in_flight {
            token_tx.send(()).expect("token channel has capacity");
        }
        let writer = scope.spawn(move || write_loop(done_rx, token_tx, sinks, progress));

        let mut dispatch_error = None;
        pool.in_place_scope(|ps| {
            let mut index = 0;
            loop {
                if token_rx.recv().is_err() {
                    break; // the writer stopped (with an error)
                }
                let next = receive(Some(&rx1)).and_then(|c1| Ok((c1, receive(rx2.as_ref())?)));
                let job = match next {
                    Err(e) => Err(e),
                    Ok((None, None)) => break,
                    Ok((Some(r1), None)) if !paired => Ok((r1, None)),
                    Ok((Some(r1), Some(r2))) if r1.len() == r2.len() => Ok((r1, Some(r2))),
                    Ok((c1, c2)) => {
                        let reads = |c: &Option<Chunk>, other: &Option<Chunk>| match (c, other) {
                            (Some(c), _) => c.first_record + c.len() as u64,
                            (None, Some(o)) => o.first_record,
                            (None, None) => 0,
                        };
                        let (n1, n2) = (reads(&c1, &c2), reads(&c2, &c1));
                        let (shorter, n) = if n1 < n2 { ("R1", n1) } else { ("R2", n2) };
                        Err(format!(
                            "R1 and R2 have different numbers of reads: {shorter} ended after {} reads",
                            group(n)
                        ))
                    }
                };
                match job {
                    Ok((r1, r2)) => {
                        let tx = done_tx.clone();
                        let job_index = index;
                        ps.spawn(move |_| {
                            let done = process(&r1, r2.as_ref()).map(|(mut outputs, stats)| {
                                for (slot, data) in outputs.iter_mut().enumerate() {
                                    if gz[slot] && !data.is_empty() {
                                        *data = gzip_member(data, level);
                                    }
                                }
                                Done {
                                    index: job_index,
                                    records: r1.len() as u64,
                                    outputs,
                                    stats,
                                }
                            });
                            let _ = tx.send(done);
                        });
                        index += 1;
                    }
                    Err(e) => {
                        dispatch_error = Some(e);
                        break;
                    }
                }
            }
        });
        drop(done_tx);
        drop(rx1);
        drop(rx2);
        let written = writer.join().expect("writer thread panicked");
        r1_thread.join().expect("reader thread panicked");
        if let Some(t) = r2_thread {
            t.join().expect("reader thread panicked");
        }
        match dispatch_error {
            Some(e) => Err(e),
            None => written,
        }
    })
}
