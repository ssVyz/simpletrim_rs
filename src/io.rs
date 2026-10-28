//! File handling: opening (optionally gzip-compressed) inputs, safe creation of outputs and
//! gzip compression of output blocks.

use std::fs::{self, File, OpenOptions};
use std::io::{BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::sync::Arc;
use std::sync::atomic::{AtomicU64, Ordering};

use flate2::Compression;
use flate2::bufread::MultiGzDecoder;
use flate2::write::GzEncoder;

use crate::Result;
use crate::fastq::FastqReader;

/// Whether a path denotes gzip data, judged by its extension.
pub fn is_gz(path: &Path) -> bool {
    path.extension()
        .is_some_and(|e| e.eq_ignore_ascii_case("gz"))
}

/// Counts the bytes read from the underlying (compressed) file, for progress reporting.
struct CountingReader<R> {
    inner: R,
    count: Arc<AtomicU64>,
}

impl<R: Read> Read for CountingReader<R> {
    fn read(&mut self, buf: &mut [u8]) -> std::io::Result<usize> {
        let n = self.inner.read(buf)?;
        self.count.fetch_add(n as u64, Ordering::Relaxed);
        Ok(n)
    }
}

/// An opened FASTQ input.
pub struct Input {
    pub reader: FastqReader,
    /// Bytes consumed from the file on disk so far.
    pub bytes_read: Arc<AtomicU64>,
    /// Size of the file on disk.
    pub file_size: u64,
}

/// Opens a FASTQ file, decompressing it if the name ends in `.gz`.
pub fn open_fastq(path: &Path, label: &str) -> Result<Input> {
    let file = File::open(path).map_err(|e| format!("cannot open {}: {e}", path.display()))?;
    let file_size = file.metadata().map(|m| m.len()).unwrap_or(0);
    let bytes_read = Arc::new(AtomicU64::new(0));
    let counted = CountingReader {
        inner: file,
        count: bytes_read.clone(),
    };
    let reader: Box<dyn std::io::BufRead + Send> = if is_gz(path) {
        let decoder = MultiGzDecoder::new(BufReader::with_capacity(256 * 1024, counted));
        Box::new(BufReader::with_capacity(1024 * 1024, decoder))
    } else {
        Box::new(BufReader::with_capacity(1024 * 1024, counted))
    };
    let label = format!("{label} ({})", path.display());
    Ok(Input {
        reader: FastqReader::new(reader, label),
        bytes_read,
        file_size,
    })
}

/// Checks that an input path exists and is a regular file.
pub fn check_input(path: &Path, what: &str) -> Result<()> {
    match fs::metadata(path) {
        Ok(m) if m.is_file() => Ok(()),
        Ok(_) => Err(format!("{what} {} is not a regular file", path.display())),
        Err(e) => Err(format!("cannot access {what} {}: {e}", path.display())),
    }
}

/// Absolute, normalised form of a path for comparing whether two names refer to the same file.
fn identity(path: &Path) -> PathBuf {
    let abs = fs::canonicalize(path)
        .or_else(|_| std::path::absolute(path))
        .unwrap_or_else(|_| path.to_path_buf());
    if cfg!(windows) {
        PathBuf::from(abs.to_string_lossy().to_lowercase())
    } else {
        abs
    }
}

/// True if both paths name the same file (best effort, also for files that do not exist yet).
pub fn same_path(a: &Path, b: &Path) -> bool {
    identity(a) == identity(b)
}

/// Validates the output paths before any work starts: all distinct, none equal to an input,
/// parent directories present, and no existing file unless `force`.
pub fn check_outputs(outputs: &[(&str, &Path)], inputs: &[&Path], force: bool) -> Result<()> {
    for (i, (what, path)) in outputs.iter().enumerate() {
        if let Some(input) = inputs.iter().find(|inp| same_path(inp, path)) {
            return Err(format!(
                "{what} {} is the same file as input {}",
                path.display(),
                input.display()
            ));
        }
        if let Some((other, _)) = outputs[..i].iter().find(|(_, p)| same_path(p, path)) {
            return Err(format!(
                "{what} and {other} both point to {}",
                path.display()
            ));
        }
        let parent = path.parent().filter(|p| !p.as_os_str().is_empty());
        if let Some(dir) = parent
            && !dir.is_dir()
        {
            return Err(format!(
                "directory {} for {what} does not exist",
                dir.display()
            ));
        }
        match fs::metadata(path) {
            Ok(m) if m.is_dir() => {
                return Err(format!("{what} {} is a directory", path.display()));
            }
            Ok(_) if !force => {
                return Err(format!(
                    "{what} {} already exists (use --force to overwrite)",
                    path.display()
                ));
            }
            _ => {}
        }
    }
    Ok(())
}

/// Output files created by this run. Files are deleted again on drop unless [`commit`] was
/// called, so a failed run does not leave partial results behind.
///
/// [`commit`]: OutputFiles::commit
pub struct OutputFiles {
    paths: Vec<PathBuf>,
    committed: bool,
}

impl OutputFiles {
    pub fn new() -> Self {
        OutputFiles {
            paths: Vec::new(),
            committed: false,
        }
    }

    /// Creates (or with `force`, truncates) an output file.
    pub fn create(&mut self, path: &Path, force: bool) -> Result<File> {
        let mut opts = OpenOptions::new();
        opts.write(true);
        if force {
            opts.create(true).truncate(true);
        } else {
            opts.create_new(true);
        }
        let file = opts.open(path).map_err(|e| {
            if e.kind() == std::io::ErrorKind::AlreadyExists {
                format!(
                    "{} already exists (use --force to overwrite)",
                    path.display()
                )
            } else {
                format!("cannot create {}: {e}", path.display())
            }
        })?;
        self.paths.push(path.to_path_buf());
        Ok(file)
    }

    /// Keeps the files.
    pub fn commit(&mut self) {
        self.committed = true;
    }
}

impl Drop for OutputFiles {
    fn drop(&mut self) {
        if !self.committed {
            for p in &self.paths {
                let _ = fs::remove_file(p);
            }
        }
    }
}

/// Compresses `data` into one complete gzip member. Members written back to back form a valid
/// (multi-member) gzip file, which lets output blocks be compressed in parallel.
pub fn gzip_member(data: &[u8], level: u32) -> Vec<u8> {
    let mut enc = GzEncoder::new(
        Vec::with_capacity(data.len() / 3 + 64),
        Compression::new(level),
    );
    enc.write_all(data).expect("writing to a Vec cannot fail");
    enc.finish().expect("writing to a Vec cannot fail")
}

/// Writes a small text file (e.g. a JSON report), honouring `force`.
pub fn write_text_file(path: &Path, content: &str, force: bool) -> Result<()> {
    let mut files = OutputFiles::new();
    let mut f = files.create(path, force)?;
    f.write_all(content.as_bytes())
        .map_err(|e| format!("cannot write {}: {e}", path.display()))?;
    files.commit();
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gz_extension() {
        assert!(is_gz(Path::new("a.fastq.gz")));
        assert!(is_gz(Path::new("a.FQ.GZ")));
        assert!(!is_gz(Path::new("a.fastq")));
    }

    #[test]
    fn gzip_members_concatenate() {
        let mut data = gzip_member(b"hello ", 4);
        data.extend(gzip_member(b"world", 4));
        let mut out = String::new();
        MultiGzDecoder::new(&data[..])
            .read_to_string(&mut out)
            .unwrap();
        assert_eq!(out, "hello world");
    }
}
