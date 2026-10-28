//! Streaming FASTQ parsing in chunks.
//!
//! A [`Chunk`] keeps the raw bytes of a batch of records in one buffer plus the positions of
//! each line, so trimming can work on slices without copying reads.

use std::io::BufRead;
use std::sync::mpsc::SyncSender;

use crate::Result;

#[derive(Clone, Copy, Debug, Default)]
struct Span {
    start: usize,
    end: usize,
}

#[derive(Clone, Copy, Debug, Default)]
struct RecordSpans {
    name: Span,
    seq: Span,
    plus: Span,
    qual: Span,
}

/// A batch of consecutive FASTQ records.
pub struct Chunk {
    /// Zero-based number of the first record in the file.
    pub first_record: u64,
    data: Vec<u8>,
    records: Vec<RecordSpans>,
}

/// A borrowed FASTQ record. `name` includes the leading `@`, `plus` the leading `+`.
#[derive(Clone, Copy, Debug)]
pub struct Record<'a> {
    pub name: &'a [u8],
    pub seq: &'a [u8],
    pub plus: &'a [u8],
    pub qual: &'a [u8],
}

impl Chunk {
    pub fn len(&self) -> usize {
        self.records.len()
    }

    pub fn get(&self, i: usize) -> Record<'_> {
        let r = &self.records[i];
        let s = |span: Span| &self.data[span.start..span.end];
        Record {
            name: s(r.name),
            seq: s(r.seq),
            plus: s(r.plus),
            qual: s(r.qual),
        }
    }

    pub fn iter(&self) -> impl Iterator<Item = Record<'_>> {
        (0..self.len()).map(|i| self.get(i))
    }
}

impl Record<'_> {
    /// Appends the record, keeping bases `[start, end)`.
    pub fn write_trimmed(&self, start: usize, end: usize, out: &mut Vec<u8>) {
        out.extend_from_slice(self.name);
        out.push(b'\n');
        out.extend_from_slice(&self.seq[start..end]);
        out.push(b'\n');
        out.extend_from_slice(self.plus);
        out.push(b'\n');
        out.extend_from_slice(&self.qual[start..end]);
        out.push(b'\n');
    }

    /// Appends the record unchanged.
    pub fn write(&self, out: &mut Vec<u8>) {
        self.write_trimmed(0, self.seq.len(), out);
    }

    /// Appends bases `[start, end)` as a record. With `piece`, `_<piece>` is appended to the
    /// read id (before any comment); with `tag`, ` <tag>` is appended to the header line.
    pub fn write_piece(
        &self,
        start: usize,
        end: usize,
        piece: Option<usize>,
        tag: Option<&str>,
        out: &mut Vec<u8>,
    ) {
        match piece {
            Some(n) => {
                let id_end = self
                    .name
                    .iter()
                    .position(|b| b.is_ascii_whitespace())
                    .unwrap_or(self.name.len());
                out.extend_from_slice(&self.name[..id_end]);
                out.extend_from_slice(format!("_{n}").as_bytes());
                out.extend_from_slice(&self.name[id_end..]);
            }
            None => out.extend_from_slice(self.name),
        }
        if let Some(tag) = tag {
            out.push(b' ');
            out.extend_from_slice(tag.as_bytes());
        }
        out.push(b'\n');
        out.extend_from_slice(&self.seq[start..end]);
        out.push(b'\n');
        out.extend_from_slice(self.plus);
        out.push(b'\n');
        out.extend_from_slice(&self.qual[start..end]);
        out.push(b'\n');
    }

    /// Read name without `@`, cut at the first whitespace.
    pub fn id(&self) -> &[u8] {
        let name = self.name.strip_prefix(b"@").unwrap_or(self.name);
        let end = name
            .iter()
            .position(|b| b.is_ascii_whitespace())
            .unwrap_or(name.len());
        &name[..end]
    }
}

/// True if two mates carry the same read name: compared up to the first whitespace, ignoring a
/// trailing `/1` or `/2`.
pub fn mate_names_match(r1: &Record, r2: &Record) -> bool {
    fn core(id: &[u8]) -> &[u8] {
        match id {
            [rest @ .., b'/', b'1' | b'2'] => rest,
            _ => id,
        }
    }
    core(r1.id()) == core(r2.id())
}

/// Reads FASTQ records from a buffered source.
pub struct FastqReader {
    inner: Box<dyn BufRead + Send>,
    label: String,
    records: u64,
    line: u64,
}

impl FastqReader {
    /// `label` names the source in error messages (e.g. `R1 (sample_R1.fastq.gz)`).
    pub fn new(inner: Box<dyn BufRead + Send>, label: String) -> Self {
        FastqReader {
            inner,
            label,
            records: 0,
            line: 0,
        }
    }

    /// Number of records returned so far.
    pub fn records_read(&self) -> u64 {
        self.records
    }

    /// Reads one line into `data`; returns its span without the line ending, or `None` at EOF.
    fn read_line(&mut self, data: &mut Vec<u8>) -> Result<Option<Span>> {
        let start = data.len();
        let n = self
            .inner
            .read_until(b'\n', data)
            .map_err(|e| format!("{}: read error: {e}", self.label))?;
        if n == 0 {
            return Ok(None);
        }
        self.line += 1;
        let mut end = data.len();
        if data[end - 1] == b'\n' {
            end -= 1;
        }
        if end > start && data[end - 1] == b'\r' {
            end -= 1;
        }
        Ok(Some(Span { start, end }))
    }

    fn required_line(&mut self, data: &mut Vec<u8>, what: &str) -> Result<Span> {
        self.read_line(data)?.ok_or_else(|| {
            format!(
                "{}: unexpected end of file in record {} (missing {what} line); the file may be truncated",
                self.label,
                self.records + 1
            )
        })
    }

    /// Reads up to `max_records` records, stopping early once the chunk holds `max_bytes`.
    /// Returns `None` at end of input.
    pub fn read_chunk(&mut self, max_records: usize, max_bytes: usize) -> Result<Option<Chunk>> {
        let mut data = Vec::with_capacity(max_bytes.min(max_records.saturating_mul(400)));
        let mut records = Vec::with_capacity(max_records.min(1 << 16));
        let first_record = self.records;
        while records.len() < max_records && data.len() < max_bytes {
            let Some(name) = self.read_line(&mut data)? else {
                break;
            };
            if name.start == name.end {
                // tolerate blank lines between records (typically at the end of the file)
                data.truncate(name.start);
                continue;
            }
            if data[name.start] != b'@' {
                return Err(format!(
                    "{}: line {}: expected a FASTQ header starting with '@' (record {})",
                    self.label,
                    self.line,
                    self.records + 1
                ));
            }
            let seq = self.required_line(&mut data, "sequence")?;
            let plus = self.required_line(&mut data, "'+' separator")?;
            if plus.start == plus.end || data[plus.start] != b'+' {
                return Err(format!(
                    "{}: line {}: expected '+' separator line (record {})",
                    self.label,
                    self.line,
                    self.records + 1
                ));
            }
            let qual = self.required_line(&mut data, "quality")?;
            if qual.end - qual.start != seq.end - seq.start {
                return Err(format!(
                    "{}: line {}: sequence and quality lengths differ ({} vs {}) in record {}",
                    self.label,
                    self.line,
                    seq.end - seq.start,
                    qual.end - qual.start,
                    self.records + 1
                ));
            }
            records.push(RecordSpans {
                name,
                seq,
                plus,
                qual,
            });
            self.records += 1;
        }
        if records.is_empty() {
            return Ok(None);
        }
        Ok(Some(Chunk {
            first_record,
            data,
            records,
        }))
    }
}

/// Reads chunks of at most `chunk_records` records (and about `chunk_bytes` bytes) and sends
/// them until the input ends, `limit` records were read (0 = no limit) or the receiver hangs
/// up. A read error is sent as the last message. Meant to run on its own thread.
pub fn feed_chunks(
    mut reader: FastqReader,
    tx: SyncSender<Result<Chunk>>,
    limit: u64,
    chunk_records: usize,
    chunk_bytes: usize,
) {
    loop {
        let wanted = if limit == 0 {
            chunk_records
        } else {
            (limit - reader.records_read()).min(chunk_records as u64) as usize
        };
        if wanted == 0 {
            return;
        }
        match reader.read_chunk(wanted, chunk_bytes) {
            Ok(Some(chunk)) => {
                if tx.send(Ok(chunk)).is_err() {
                    return;
                }
            }
            Ok(None) => return,
            Err(e) => {
                let _ = tx.send(Err(e));
                return;
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    fn reader(text: &str) -> FastqReader {
        FastqReader::new(
            Box::new(Cursor::new(text.as_bytes().to_vec())),
            "test".into(),
        )
    }

    #[test]
    fn parses_records_in_chunks() {
        let mut r = reader("@r1 a\nACGT\n+\nIIII\r\n@r2\nGG\n+r2\nII\n\n@r3\nT\n+\nI");
        let c = r.read_chunk(2, usize::MAX).unwrap().unwrap();
        assert_eq!(c.len(), 2);
        assert_eq!(c.first_record, 0);
        let rec = c.get(0);
        assert_eq!(rec.name, b"@r1 a");
        assert_eq!(rec.qual, b"IIII");
        assert_eq!(rec.id(), b"r1");
        assert_eq!(c.get(1).plus, b"+r2");
        let c = r.read_chunk(2, usize::MAX).unwrap().unwrap();
        assert_eq!(c.len(), 1);
        assert_eq!(c.first_record, 2);
        assert_eq!(c.get(0).seq, b"T");
        assert!(r.read_chunk(2, usize::MAX).unwrap().is_none());
        assert_eq!(r.records_read(), 3);
    }

    #[test]
    fn writes_pieces_with_suffix_and_tag() {
        let mut r = reader("@id\tqs:f:12.5 x\nACGTACGT\n+\nIIIIHHHH\n");
        let c = r.read_chunk(1, usize::MAX).unwrap().unwrap();
        let mut out = Vec::new();
        c.get(0).write_piece(2, 6, Some(2), None, &mut out);
        c.get(0)
            .write_piece(0, 2, None, Some("failed_too_short"), &mut out);
        assert_eq!(
            String::from_utf8(out).unwrap(),
            "@id_2\tqs:f:12.5 x\nGTAC\n+\nIIHH\n@id\tqs:f:12.5 x failed_too_short\nAC\n+\nII\n"
        );
    }

    #[test]
    fn byte_limited_chunks() {
        let mut r = reader("@a\nACGT\n+\nIIII\n@b\nACGT\n+\nIIII\n@c\nA\n+\nI\n");
        // the chunk ends once it holds at least 20 bytes
        assert_eq!(r.read_chunk(100, 20).unwrap().unwrap().len(), 2);
        assert_eq!(r.read_chunk(100, 20).unwrap().unwrap().len(), 1);
    }

    #[test]
    fn rejects_malformed_input() {
        assert!(
            reader("r1\nACGT\n+\nIIII\n")
                .read_chunk(10, usize::MAX)
                .is_err()
        );
        assert!(
            reader("@r1\nACGT\n-\nIIII\n")
                .read_chunk(10, usize::MAX)
                .is_err()
        );
        assert!(
            reader("@r1\nACGT\n+\nIII\n")
                .read_chunk(10, usize::MAX)
                .is_err()
        );
        assert!(reader("@r1\nACGT\n+\n").read_chunk(10, usize::MAX).is_err());
    }

    #[test]
    fn writes_trimmed_records() {
        let mut r = reader("@r1\nACGTAC\n+\nABCDEF\n");
        let c = r.read_chunk(10, usize::MAX).unwrap().unwrap();
        let mut out = Vec::new();
        c.get(0).write_trimmed(1, 4, &mut out);
        assert_eq!(out, b"@r1\nCGT\n+\nBCD\n");
    }

    #[test]
    fn mate_names() {
        let mut r = reader(
            "@a/1\nA\n+\nI\n@a/2\nA\n+\nI\n@x 1:N:0\nA\n+\nI\n@x 2:N:0\nA\n+\nI\n@y\nA\n+\nI\n",
        );
        let c = r.read_chunk(10, usize::MAX).unwrap().unwrap();
        assert!(mate_names_match(&c.get(0), &c.get(1)));
        assert!(mate_names_match(&c.get(2), &c.get(3)));
        assert!(!mate_names_match(&c.get(0), &c.get(4)));
    }
}
