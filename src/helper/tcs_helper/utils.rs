use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::{Result as IoResult, Write};
use std::ops::Range;

use bio::io::fastq::Record;
use chrono::Local;

pub fn log_line(writer: &mut BufWriter<File>, message: &str) -> IoResult<()> {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S");
    writeln!(writer, "[{}] {}", now, message)?;
    writer.flush()?;
    Ok(())
}

pub fn diff_positions(a: &str, b: &str) -> Vec<usize> {
    a.chars()
        .zip(b.chars())
        .enumerate()
        .filter_map(|(i, (ca, cb))| if ca != cb { Some(i) } else { None })
        .collect()
}

static IUPAC_TUPLES: &'static [(char, &'static [char])] = &[
    ('A', &['A']),
    ('C', &['C']),
    ('G', &['G']),
    ('T', &['T']),
    ('R', &['A', 'G']),
    ('Y', &['C', 'T']),
    ('S', &['G', 'C']),
    ('W', &['A', 'T']),
    ('K', &['G', 'T']),
    ('M', &['A', 'C']),
    ('B', &['C', 'G', 'T']),
    ('D', &['A', 'G', 'T']),
    ('H', &['A', 'C', 'T']),
    ('V', &['A', 'C', 'G']),
    ('N', &['A', 'C', 'G', 'T']),
];

pub fn get_iupac_bases(c: char) -> Option<&'static [char]> {
    let c = c.to_ascii_uppercase();
    IUPAC_TUPLES
        .iter()
        .find_map(|&(base, bases)| if base == c { Some(bases) } else { None })
}

pub fn iupac_matches(a: char, b: char) -> bool {
    if a == b {
        return true;
    }

    if let (Some(a_bases), Some(b_bases)) = (get_iupac_bases(a), get_iupac_bases(b)) {
        a_bases.iter().any(|base| b_bases.contains(base))
    } else {
        false
    }
}

pub fn diff_by_iupac(a: &str, b: &str) -> Vec<usize> {
    a.chars()
        .zip(b.chars())
        .enumerate()
        .filter_map(
            |(i, (ca, cb))| {
                if iupac_matches(ca, cb) { None } else { Some(i) }
            },
        )
        .collect()
}

pub fn diff_byte_equal_length(a: &[u8], b: &[u8]) -> Vec<usize> {
    (0..a.len()).filter(|&i| a[i] != b[i]).collect()
}

pub trait FastqRecordTrimExt {
    /// Trims the read and quality to the specified length.
    fn get_range(&self, length: Range<usize>) -> Result<Record, Box<dyn Error + Send + Sync>>;
}

impl FastqRecordTrimExt for Record {
    fn get_range(&self, length: Range<usize>) -> Result<Record, Box<dyn Error + Send + Sync>> {
        if length.start >= self.seq().len() || length.end > self.seq().len() {
            return Err("Invalid range for trimming: start or end exceeds sequence length.".into());
        } else {
            let seq = &self.seq()[length.clone()];
            let qual = &self.qual()[length];

            Ok(Record::with_attrs(self.id(), self.desc(), seq, qual))
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_diff_positions() {
        let a = "ACGT";
        let b = "AGCT";
        let diff = diff_positions(a, b);
        assert_eq!(diff, vec![1, 2]);
    }

    #[test]
    fn test_iupac_matches() {
        assert!(iupac_matches('A', 'A'));
        assert!(iupac_matches('R', 'A'));
        assert!(!iupac_matches('A', 'C'));
    }

    #[test]
    fn test_diff_by_iupac() {
        let a = "ACGTRC";
        let b = "AGCTGW";
        let diff = diff_by_iupac(a, b);
        assert_eq!(diff, vec![1, 2, 5]);
    }
}
