use std::error::Error;
use std::fs::File;
use std::io::BufWriter;
use std::io::{Result as IoResult, Write};
use std::ops::Range;

use bio::alphabets::dna;
use bio::io::fasta;
use bio::io::fastq::{self, Record};
use chrono::Local;
use virust_locator::prelude::*;

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

pub fn trim_sequence_from_locator(
    locator: &Locator,
    start: usize,
    end: usize,
) -> Result<(Vec<u8>, Range<usize>), Box<dyn Error + Send + Sync>> {
    let query_aligned = locator.query_aligned_string.as_bytes();
    let ref_aligned = locator.ref_aligned_string.as_bytes();

    let mut l1 = locator.ref_start;
    let mut l2 = locator.ref_end;

    let mut g1 = 0;
    let mut g2 = 0;

    for c in ref_aligned.iter() {
        if l1 == start {
            break;
        }
        g1 += 1;
        if *c != b'-' {
            // '-' is 45 in ASCII
            l1 += 1;
        }
    }

    for c in ref_aligned.iter().rev() {
        if l2 == end {
            break;
        }
        g2 += 1;
        if *c != b'-' {
            l2 -= 1;
        }
    }

    // TODO: implement error handling for cases where g1 or g2 exceed the length of query_aligned
    if g1 >= query_aligned.len() || g2 >= query_aligned.len() || (g1 + g2) > query_aligned.len() {
        return Err("Gaps exceed the length of the aligned query sequence.".into());
    }

    let trimmed_seq = query_aligned[g1..query_aligned.len() - g2]
        .iter()
        .cloned()
        .filter(|&c| c != b'-') // Remove gaps
        .collect();

    let prefix = query_aligned[..g1].iter().filter(|&&c| c != b'-').count();
    let suffix = query_aligned[query_aligned.len() - g2..]
        .iter()
        .filter(|&&c| c != b'-')
        .count();

    let trimmed_range = prefix..(query_aligned.len() - suffix);
    Ok((trimmed_seq, trimmed_range))
}

pub fn reverse_complement(record: &Record) -> Record {
    let seq = record
        .seq()
        .iter()
        .rev()
        .map(|&c| dna::complement(c))
        .collect::<Vec<u8>>();

    let qual = record.qual().iter().rev().cloned().collect::<Vec<u8>>();

    Record::with_attrs(record.id(), record.desc(), &seq, &qual)
}

pub fn fastq_to_fasta_record(fq: &fastq::Record) -> fasta::Record {
    fasta::Record::with_attrs(fq.id(), fq.desc(), fq.seq())
}

pub static CLI_ANIMATION_TICK_STRINGS: [&str; 26] = [
    "🐱  🐭          ",
    " 🐱  🐭         ",
    "  🐱  🐭        ",
    "   🐱  🐭       ",
    "    🐱  🐭      ",
    "     🐱  🐭     ",
    "      🐱  🐭    ",
    "       🐱  🐭   ",
    "        🐱  🐭  ",
    "         🐱  🐭 ",
    "          🐱  🐭",
    "           🐱🐭 ",
    "           🐭 🐱",
    "         🐭  🐱 ",
    "        🐭  🐱  ",
    "       🐭  🐱   ",
    "      🐭  🐱    ",
    "     🐭  🐱     ",
    "    🐭  🐱      ",
    "   🐭  🐱       ",
    "  🐭  🐱        ",
    " 🐭  🐱         ",
    "🐭  🐱          ",
    " 🐭🐱           ",
    "🐱🐭            ",
    "    🐱 ❤️  🐭    ",
];

/// This constant defines the threshold for low abundance in raw reads.
/// It is set to 0.0005, which means that if the abundance of reads from one Region is less than 0.05% of the total raw reads, it will be considered low abundance.
/// Potentially from cross-contamination in library preparation or sequencing.
pub const LOW_ABUNDANCE_THRESHOLD_FOR_RAW_READS: f64 = 0.0005; // 0.05% of raw reads

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

    #[test]

    fn test_trim_sequence_from_locator() {
        let seq = "ATCCTTTAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCTATTGAGACTGTACCAGTAA";
        let pr = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT".to_string();
        let args = Args {
            query: vec![seq.to_string()],
            reference: "HXB2".to_string(),
            type_query: "nt".to_string(),
            algorithm: 1,
        };

        let locator = Locator::build(&args).unwrap().pop().unwrap().unwrap();

        let trimmed = trim_sequence_from_locator(&locator, 2253, 2549).unwrap();

        let trimmed_str = String::from_utf8_lossy(&trimmed.0);

        let trimmed_range = trimmed.1;

        dbg!(&trimmed_str);
        dbg!(&trimmed_range);
        assert_eq!(trimmed_str.to_string(), pr);
    }

    #[test]
    fn test_trim_sequence_from_locator_with_gaps() {
        let seq = "ATCCAACTTCCCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTTCCCATTAGCCCGACTATTGAGACTGTACCAGTAA";
        let pr = "CCTCAGGTCACTCTTTGGCAACGACCCCTCGTCACAATAAAGATAGGGGGGCAACTAAAGGAAGCTCTATTAGATACAGGAGCAGATGATACAGTATTAGAAGAAATGAGTTTGCCAGGAAGATGGAAACCAAAAATGATAGGGGGAATTGGAGGTTTTATCAAAGTAAGACAGTATGATCAGATACTCATAGAAATCTGTGGACATAAAGCTATAGGTACAGTATTAGTAGGACCTACACCTGTCAACATAATTGGAAGAAATCTGTTGACTCAGATTGGTTGCACTTTAAATTTT".to_string();
        let args = Args {
            query: vec![seq.to_string()],
            reference: "HXB2".to_string(),
            type_query: "nt".to_string(),
            algorithm: 1,
        };

        let locator = Locator::build(&args).unwrap().pop().unwrap().unwrap();

        let trimmed = trim_sequence_from_locator(&locator, 2253, 2549).unwrap();

        let trimmed_str = String::from_utf8_lossy(&trimmed.0);

        let trimmed_range = trimmed.1;

        assert_eq!(trimmed_str.to_string(), pr);
        assert_eq!(trimmed_range, 10..307);
    }

    #[test]
    fn test_reverse_complement() {
        let record = Record::with_attrs("test", None, b"ATCG", b"1234");

        let rev_comp = reverse_complement(&record);

        assert_eq!(rev_comp.seq(), b"CGAT");
        assert_eq!(rev_comp.qual(), b"4321");
        assert_eq!(rev_comp.id(), "test");
        assert_eq!(rev_comp.desc(), None);
    }
}
