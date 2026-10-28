//! Guessing the sequencing platform from FASTQ read headers.

use serde::Serialize;

#[derive(Clone, Copy, Debug, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Platform {
    Illumina,
    ElementAviti,
    Mgi,
    OxfordNanopore,
    PacBio,
    IonTorrent,
    Unknown,
}

impl Platform {
    pub fn label(self) -> &'static str {
        match self {
            Platform::Illumina => "Illumina",
            Platform::ElementAviti => "Element Biosciences AVITI",
            Platform::Mgi => "MGI / DNBSEQ",
            Platform::OxfordNanopore => "Oxford Nanopore",
            Platform::PacBio => "PacBio",
            Platform::IonTorrent => "Ion Torrent",
            Platform::Unknown => "unknown",
        }
    }

    pub fn is_long_read(self) -> bool {
        matches!(self, Platform::OxfordNanopore | Platform::PacBio)
    }
}

/// What a read header tells about its origin.
#[derive(Clone, Debug, PartialEq, Eq, Serialize)]
pub struct HeaderInfo {
    pub platform: Platform,
    /// Header layout, e.g. `casava_1.8`, `illumina_legacy`, `mgi`, `ont`, `sra`.
    pub format: &'static str,
    pub instrument_id: Option<String>,
    pub instrument_model: Option<&'static str>,
    /// Mate number stated in the header (`1:N:...` or `/1`), if any.
    pub read_number: Option<u8>,
    /// Instrument uses two-colour chemistry, where no signal is called as G (poly-G tails).
    pub two_color: bool,
}

/// Instrument id prefixes of two-colour instruments (from fastp, see fastp PR #508). fastp
/// matches a bare `@A` for the NovaSeq 6000; here `A` must be followed by a digit so that
/// e.g. Element AVITI (`AV...`) is not mistaken for a NovaSeq.
const TWO_COLOR_PREFIXES: &[&str] = &["FS", "MN", "SH", "NS", "NB", "NDX", "VL", "VH", "NA", "LH"];

/// True if the read name comes from a two-colour instrument (fastp's polyG auto-enable rule).
pub fn is_two_color(name: &[u8]) -> bool {
    let id = name.strip_prefix(b"@").unwrap_or(name);
    if TWO_COLOR_PREFIXES
        .iter()
        .any(|p| id.starts_with(p.as_bytes()))
    {
        return true;
    }
    matches!(id, [b'A', d, ..] if d.is_ascii_digit())
}

fn starts_with_then_digit(s: &str, prefix: &str) -> bool {
    s.strip_prefix(prefix)
        .and_then(|rest| rest.bytes().next())
        .is_some_and(|b| b.is_ascii_digit())
}

/// Instrument model from an Illumina-style instrument id (best effort).
fn instrument_model(id: &str) -> (Platform, Option<&'static str>) {
    const TABLE: &[(&str, &str)] = &[
        ("LH", "NovaSeq X"),
        ("VH", "NextSeq 1000/2000"),
        ("VL", "NextSeq 1000/2000"),
        ("NDX", "NextSeq 550Dx"),
        ("NB", "NextSeq 500/550"),
        ("NS", "NextSeq 500/550"),
        ("NA", "NovaSeq"),
        ("MN", "MiniSeq"),
        ("SH", "MiniSeq"),
        ("FS", "iSeq 100"),
        ("A", "NovaSeq 6000"),
        ("M", "MiSeq"),
        ("K", "HiSeq 3000/4000"),
        ("J", "HiSeq 3000"),
        ("E", "HiSeq X"),
        ("D", "HiSeq 2500"),
        ("C", "HiSeq 1500"),
    ];
    if starts_with_then_digit(id, "AV") {
        return (Platform::ElementAviti, Some("AVITI"));
    }
    if id.starts_with("HWI-M") {
        return (Platform::Illumina, Some("MiSeq"));
    }
    if id.starts_with("HWI-ST") || id.starts_with("HWI-D") {
        return (Platform::Illumina, Some("HiSeq 2000/2500"));
    }
    if id.starts_with("HWUSI") {
        return (Platform::Illumina, Some("Genome Analyzer IIx"));
    }
    let model = TABLE
        .iter()
        .find(|(p, _)| starts_with_then_digit(id, p))
        .map(|(_, m)| *m);
    (Platform::Illumina, model)
}

fn is_uuid(s: &str) -> bool {
    let b = s.as_bytes();
    b.len() == 36
        && b.iter().enumerate().all(|(i, &c)| match i {
            8 | 13 | 18 | 23 => c == b'-',
            _ => c.is_ascii_hexdigit(),
        })
}

fn all_digits(s: &str) -> bool {
    !s.is_empty() && s.bytes().all(|b| b.is_ascii_digit())
}

/// `V350012345L1C001R0010000001`-style DNBSEQ read ids.
fn is_mgi_id(id: &str) -> bool {
    let b = id.as_bytes();
    let letters = b.iter().take_while(|c| c.is_ascii_uppercase()).count();
    if !(1..=2).contains(&letters) {
        return false;
    }
    let rest = &id[letters..];
    let (Some(l), Some(c), Some(r)) = (rest.find('L'), rest.find('C'), rest.find('R')) else {
        return false;
    };
    l < c
        && c < r
        && all_digits(&rest[..l])
        && all_digits(&rest[l + 1..c])
        && all_digits(&rest[c + 1..r])
        && {
            let tail = rest[r + 1..].trim_end_matches(['/', '1', '2']);
            !tail.is_empty() && tail.bytes().all(|b| b.is_ascii_digit() || b == b'_')
        }
}

fn slash_mate(id: &str) -> Option<u8> {
    if id.ends_with("/1") {
        Some(1)
    } else if id.ends_with("/2") {
        Some(2)
    } else {
        None
    }
}

/// Parses one FASTQ header line (with or without the leading `@`).
pub fn parse_header(name: &[u8]) -> HeaderInfo {
    let text = String::from_utf8_lossy(name.strip_prefix(b"@").unwrap_or(name));
    let mut tokens = text.split_ascii_whitespace();
    let first = tokens.next().unwrap_or("");
    let second = tokens.next();
    let two_color = is_two_color(name);
    let info = |platform, format, instrument_id: Option<&str>, model, read_number| HeaderInfo {
        platform,
        format,
        instrument_id: instrument_id.map(str::to_string),
        instrument_model: model,
        read_number,
        two_color,
    };

    if is_uuid(first)
        || text.contains("runid=")
        || text.contains("start_time=")
        || text.contains("parent_read_id=")
    {
        return info(Platform::OxfordNanopore, "ont", None, None, None);
    }

    let fields: Vec<&str> = first.split(':').collect();
    if (fields.len() == 7 || fields.len() == 8) && fields[3..7].iter().all(|f| all_digits(f)) {
        let read_number = second
            .and_then(|s| s.split(':').next())
            .and_then(|n| n.parse::<u8>().ok())
            .filter(|n| (1..=3).contains(n));
        let (platform, model) = instrument_model(fields[0]);
        return info(platform, "casava_1.8", Some(fields[0]), model, read_number);
    }

    if ["SRR", "ERR", "DRR"].iter().any(|p| first.starts_with(p)) {
        return info(Platform::Unknown, "sra", None, None, slash_mate(first));
    }

    if first.starts_with('m') && first.contains('/') && first.contains('_') {
        let run = first.split('/').next().unwrap_or("");
        if run[1..]
            .split('_')
            .next()
            .is_some_and(|s| s.chars().any(|c| c.is_ascii_digit()))
        {
            return info(Platform::PacBio, "pacbio", None, None, None);
        }
    }

    if is_mgi_id(first) {
        return info(Platform::Mgi, "mgi", None, None, slash_mate(first));
    }

    let core = first.split(['#', '/']).next().unwrap_or("");
    let legacy: Vec<&str> = core.split(':').collect();
    if legacy.len() == 5 && legacy[1..].iter().all(|f| all_digits(f)) {
        let (platform, model) = instrument_model(legacy[0]);
        return info(
            platform,
            "illumina_legacy",
            Some(legacy[0]),
            model,
            slash_mate(first),
        );
    }

    if fields.len() == 3 && fields[1..].iter().all(|f| all_digits(f)) {
        return info(Platform::IonTorrent, "ion_torrent", None, None, None);
    }

    info(Platform::Unknown, "unknown", None, None, slash_mate(first))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn element_aviti_is_not_two_color() {
        let h =
            parse_header(b"@AV243904:0107_A:2530453498:2:10901:0188:0019 1:N:0:GCCACAGG+CATGCCAT");
        assert_eq!(h.platform, Platform::ElementAviti);
        assert_eq!(h.format, "casava_1.8");
        assert_eq!(h.read_number, Some(1));
        assert!(!h.two_color);
    }

    #[test]
    fn novaseq_is_two_color() {
        let h = parse_header(b"@A00123:8:H5KNVDSXX:1:1101:1000:1000 2:N:0:ACGT");
        assert_eq!(h.platform, Platform::Illumina);
        assert_eq!(h.instrument_model, Some("NovaSeq 6000"));
        assert_eq!(h.read_number, Some(2));
        assert!(h.two_color);
        assert!(is_two_color(b"@LH00123:1:X:1:1:1:1"));
        assert!(!is_two_color(b"@M01234:1:X:1:1:1:1"));
    }

    #[test]
    fn miseq_legacy_and_others() {
        assert_eq!(
            parse_header(b"@M01234:55:000000000-A1B2C:1:1101:15589:1331 1:N:0:1").instrument_model,
            Some("MiSeq")
        );
        let legacy = parse_header(b"@HWUSI-EAS100R:6:73:941:1973#0/1");
        assert_eq!(legacy.format, "illumina_legacy");
        assert_eq!(legacy.read_number, Some(1));
        assert_eq!(
            parse_header(b"@V350012345L1C001R0010000001/1").platform,
            Platform::Mgi
        );
        assert_eq!(
            parse_header(b"@0a1b2c3d-1234-4abc-9def-0123456789ab runid=abc read=1 ch=5").platform,
            Platform::OxfordNanopore
        );
        assert_eq!(
            parse_header(b"@m64011_190830_220126/1/ccs").platform,
            Platform::PacBio
        );
        assert_eq!(parse_header(b"@SRR1234567.1 1 length=100").format, "sra");
        assert_eq!(parse_header(b"@whatever").platform, Platform::Unknown);
    }
}
