use std::error::Error;

use bio::io::fasta;
use bio::io::fastq;
use getset::{Getters, Setters};

const MIN_OVERLAP: usize = 10; // minimum overlap length, can be adjusted
const ERROR_RATE_FOR_ENDJOINING: f64 = 0.02; // allowed error rate, can be adjusted

/// Strategy for joining two ends of sequences.
/// This enum defines how the end joining should be performed based on the overlap information.
/// - `Simple`: No overlap check, just concatenate the sequences.
/// - `Overlap(usize)`: Join with a known overlap length.
/// - `UnknownOverlap`: Attempt to find the best overlap automatically.
/// The `Overlap` variant allows specifying a fixed overlap length, while `UnknownOverlap` will
/// try to determine the best overlap based on the sequences provided.
/// The `Simple` variant is useful when the sequences are known to be non-overlapping or when
/// the overlap is not relevant for the joining process.
/// Please note that the `Overlap` variant requires the user to know the exact overlap length
/// between the two sequences, and assumes that the sequences overlap at the end of r1 and the start of r2.
#[derive(Debug, Clone)]
pub enum EndJoiningStrategy {
    // simple joining, no overlap check
    Simple,
    // joining with known overlap,
    Overlap(usize),
    // unknown overlap, will try to find the best overlap
    UnknownOverlap,
}

/// Input for the end joining process.
/// This enum allows for different types of sequence records to be used as input.
/// - `Fasta`: A tuple containing two slices of fasta records, r1 and r2.
/// - `Fastq`: A tuple containing two slices of fastq records, r1 and r2.
/// The `Fasta` variant is used when the sequences are in FASTA format, while the `Fastq` variant
/// is used when the sequences are in FASTQ format.
/// The choice of input type allows the end joining function to handle both types of sequence data,
/// which is useful for applications that may work with either format.
#[derive(Debug, Clone)]
pub enum EndJoiningInput<'a> {
    Fasta((&'a fasta::Record, &'a fasta::Record)),
    Fastq((&'a fastq::Record, &'a fastq::Record)),
}

impl EndJoiningInput<'_> {
    /// Validates the records in the end joining input.
    pub fn validate_records(&self) -> Result<(), Box<dyn Error + Send + Sync>> {
        match self {
            EndJoiningInput::Fasta((r1, r2)) => {
                if r1.is_empty() || r2.is_empty() {
                    return Err("Fasta records cannot be empty".into());
                }
                Ok(())
            }
            EndJoiningInput::Fastq((r1, r2)) => {
                if r1.is_empty() || r2.is_empty() {
                    return Err("Fastq records cannot be empty".into());
                }
                Ok(())
            }
        }
    }
}
/// Result of the end joining process.
/// This struct contains the joined sequence and optionally the quality scores.
/// - `seq`: The joined sequence as a vector of bytes.
/// - `quality`: An optional vector of quality scores corresponding to the joined sequence.
/// The `EndJoiningResult` struct is used to represent the outcome of the end joining operation.
/// It provides the joined sequence and, if available, the quality scores.
#[derive(Debug, Clone, Getters, Setters)]
pub struct EndJoiningResult {
    #[getset(get = "pub")]
    seq: Vec<u8>,
    #[getset(get = "pub")]
    quality: Option<Vec<u8>>,
}

impl EndJoiningResult {
    /// Creates a new `EndJoiningResult` with an empty sequence and no quality scores.
    pub fn new() -> Self {
        EndJoiningResult {
            seq: Vec::new(),
            quality: None,
        }
    }
}

/// Result of an overlap between two sequences.
/// This struct contains the offset, overlap length, and number of mismatches.
/// - `offset`: The offset of the second sequence relative to the r1, the first sequence.
/// - `overlap_len`: The length of the overlap between the two sequences.
/// - `mismatches`: The number of mismatches in the overlapping region.
/// The `OverlapResult` struct is used to represent the result of an overlap analysis between two sequences.
/// It provides the offset of the second sequence relative to the first, the length of the overlap,
/// and the number of mismatches found in the overlapping region.
/// This information is crucial for understanding how the two sequences align and where they can be joined.
/// The offset is interpreted as follows:
/// # Offset interpretation:
/// - `offset == 0`: r2 aligns to the very first position of r1.
/// - `offset > 0`: r2 is shifted to the right, aligning with a later position in r1.
/// - `offset < 0`: r2 hangs off the left end of r1.
/// - `offset == r1.len()`: r2 is completely after r1; no overlap.
/// In all cases, only offsets that yield an overlap of at least `min_overlap` are considered.
/// If there is no acceptable overlap, returns None.
/// # Example: Overlap offset illustration, positive offset
/// ```text
/// pos: 0 1 2 3 4 5 6 7 8 9
/// r1:  A C G T A C G T
/// r2:        T A C G T T G
/// pos:       0 1 2 3 4 5 6
///
/// An offset of 3 means r2 starts at position 3 of r1:
/// ```
/// Overlap length: 5 (from `r1[3..8]` vs `r2[0..5]`)
///
/// # Example: Negative offset
/// ```text
/// pos:       0 1 2 3 4 5 6 7
/// r1:        A C G T A C G A
/// r2:  A G T A C G T A
/// pos: 0 1 2 3 4 5 6 7
///
/// An offset of -3 means r1 starts at position 3 of r2:
/// ```
/// Overlap length: 5 (from `r1[0..5]` vs `r2[3..8]`)
/// # Example: No overlap, aka overlap_len = 0
/// ``` text
/// pos: 0 1 2 3 4 5 6 7
/// r1:  A C G T A C G T
/// r2:                   T A C G T T G
/// pos:                  0 1 2 3 4 5 6
/// ```
/// In this case, the offset would be 8 (`r1.len()`), and the overlap length would be 0.
/// This indicates that r2 starts after r1 ends, with no overlap.
#[derive(Debug, Clone, PartialEq, Getters, Setters)]
pub struct OverlapResult {
    #[getset(get = "pub")]
    offset: isize,
    #[getset(get = "pub")]
    overlap_len: usize,
    #[getset(get = "pub")]
    mismatches: usize,
}

impl OverlapResult {
    /// Create a new OverlapResult with default values.
    pub fn new() -> Self {
        OverlapResult {
            offset: 0,
            overlap_len: 0,
            mismatches: 0,
        }
    }

    /// Create an OverlapResult from two sequences with a known overlap length.
    /// This assumes that the sequences overlap at the end of r1 and the start of r2.
    /// The offset is calculated as the length of r1 minus the overlap length.
    pub fn from_simple_overlap(len1: usize, len2: usize, overlap_len: usize) -> Self {
        // the overlap length should not exceed the length of r1 or r2
        // do not panic, just return the maximum possible overlap
        // this also ensure the offset is always non-negative
        let overlap_len = [overlap_len, len1, len2].into_iter().min().unwrap();

        let offset = len1 - overlap_len;

        OverlapResult {
            offset: offset as isize,
            overlap_len,
            mismatches: 0, // no mismatches in this case
        }
    }
}

/// End joining function that takes two sequences and joins them based on the specified strategy.
/// This function handles both FASTA and FASTQ formats, allowing for flexible input.
/// It uses the `EndJoiningInput` enum to determine the input type and the `EndJoiningStrategy` enum
/// to determine how to join the sequences.
/// The function returns an `EndJoiningResult` containing the joined sequence and optionally the quality scores.
/// # Arguments
/// - `input`: An `EndJoiningInput` enum that specifies the input type (Fasta or Fastq).
/// - `strategy`: An `EndJoiningStrategy` enum that specifies how to join the sequences (Simple, Overlap, or UnknownOverlap).
/// # Returns
/// - `Result<EndJoiningResult, Box<dyn Error + Send + Sync>>`: The result of the end joining operation.
///   - On success, it returns an `EndJoiningResult` containing the joined sequence and quality scores.
///   - On failure, it returns an error wrapped in a `Box<dyn Error + Send + Sync>`.
/// # Example
/// ```ignore
/// let input = EndJoiningInput::Fasta((fasta_records1, fasta_records2));
/// let strategy = EndJoiningStrategy::UnknownOverlap;
/// let result = end_joining(input, strategy);
/// match result {
///     Ok(joined_result) => {
///         println!("Joined sequence: {:?}", joined_result.seq);
///         if let Some(quality) = joined_result.quality {
///             println!("Quality scores: {:?}", quality);
///         } else {
///             println!("No quality scores available.");
///         }
///     }
///     Err(e) => {
///         eprintln!("Error during end joining: {}", e);
///     }
/// }
/// ```
pub fn end_joining(
    input: EndJoiningInput,
    strategy: &EndJoiningStrategy,
) -> Result<EndJoiningResult, Box<dyn Error + Send + Sync>> {
    // Validate the input records
    input.validate_records()?;
    let (r1, r2, q1, q2): (Vec<u8>, Vec<u8>, Option<Vec<u8>>, Option<Vec<u8>>) = match input {
        EndJoiningInput::Fasta((r1, r2)) => {
            let seq1 = r1.seq().to_owned();
            let seq2 = r2.seq().to_owned();
            (seq1, seq2, None, None) // Fasta does not have quality scores
        }
        EndJoiningInput::Fastq((r1, r2)) => {
            let seq1 = r1.seq().to_owned();
            let seq2 = r2.seq().to_owned();
            let qual1 = r1.qual().to_owned();
            let qual2 = r2.qual().to_owned();
            (seq1, seq2, Some(qual1), Some(qual2))
        }
    };

    let overlap = match strategy {
        EndJoiningStrategy::Simple => {
            // this is equivalent to zero overlap.
            OverlapResult::from_simple_overlap(r1.len(), r2.len(), 0)
        }
        EndJoiningStrategy::Overlap(overlap_len) => {
            // use the provided overlap length
            OverlapResult::from_simple_overlap(r1.len(), r2.len(), *overlap_len)
        }
        EndJoiningStrategy::UnknownOverlap => {
            // find the best overlap
            find_best_overlap(&r1, &r2, MIN_OVERLAP, ERROR_RATE_FOR_ENDJOINING)
        }
    };

    // Join the sequences based on the overlap result
    let end_joining_result = join_with_overlap(&r1, q1.as_deref(), &r2, q2.as_deref(), overlap);

    Ok(end_joining_result)
}

/// Finds the best overlap between two sequences based on a minimum overlap length and an error rate.
/// This function iterates through possible offsets to find the best overlap region.
/// It calculates the number of mismatches in the overlapping region and compares it to the allowed error rate.
/// If a better overlap is found, it updates the best overlap result.
/// It restrict the left overhang because r2 is supposed to be the downstream end of r1. If there is a long left overhang, it is likely a very short insert for a mistake.
/// The function returns an `OverlapResult` containing the best overlap found.
/// # Arguments
/// - `r1`: A slice of bytes representing the first sequence.
/// - `r2`: A slice of bytes representing the second sequence.
/// - `min_overlap`: The minimum length of the overlap required.
/// - `error_rate`: The allowed error rate as a fraction of the overlap length.
/// # Returns
/// - `OverlapResult`: The best overlap found between the two sequences.
///   - If no overlap is found, it returns an `OverlapResult` with `offset` set to the length of `r1`, `overlap_len` set to 0, and `mismatches` set to 0.
pub fn find_best_overlap(
    r1: &[u8],
    r2: &[u8],
    min_overlap: usize,
    error_rate: f64,
) -> OverlapResult {
    let len1 = r1.len() as isize;
    let len2 = r2.len() as isize;
    let mut best: Option<OverlapResult> = None;

    // Offsets: how much r2 is shifted relative to r1
    let half_len2 = (len2 / 2) as isize;
    let raw_min_offset = -(len2 - min_overlap as isize);
    let min_offset = raw_min_offset.max(-half_len2); // restrict left overhang
    let max_offset = len1 - min_overlap as isize; // right overhang can go up to the end of r1 minus min_overlap

    for offset in min_offset..=max_offset {
        // Determine the overlap region in r1 and r2
        let start1 = offset.max(0) as usize; // r1 starts
        let start2 = (-offset).max(0) as usize; // r2 starts
        let end1 = len1.min(offset + len2) as usize; // r1 ends
        let overlap = end1.saturating_sub(start1); // length of the overlap

        if overlap >= min_overlap
            && (start2 + overlap) <= len2 as usize
            && (start1 + overlap) <= len1 as usize
        {
            let mismatches = r1[start1..end1]
                .iter()
                .zip(r2[start2..start2 + overlap].iter())
                .filter(|(a, b)| a != b)
                .count();
            if (mismatches as f64) <= (overlap as f64 * error_rate) {
                // favor longer overlaps;
                let is_better = match &best {
                    None => true,
                    Some(best_overlap) => {
                        overlap > best_overlap.overlap_len
                            || (overlap == best_overlap.overlap_len
                                && mismatches < best_overlap.mismatches)
                    }
                };
                if is_better {
                    best = Some(OverlapResult {
                        offset,
                        overlap_len: overlap,
                        mismatches,
                    });
                }
            }
        }
    }

    if best.is_some() {
        best.unwrap()
    } else {
        OverlapResult {
            offset: len1,   // no overlap found, set offset to len1
            overlap_len: 0, // overlap_len == 0 means there is no overlap.
            mismatches: 0,
        }
    }
}

/// Joins two sequences with an overlap based on the provided `OverlapResult`.
/// This function handles the case where there is no overlap by simply concatenating the sequences.
/// If there is an overlap, it creates a consensus sequence based on the overlapping region.
/// It uses the quality scores from both sequences to determine the consensus base in the overlap region if available,
/// base with higher quality are returned as the consensus base.
/// It also builds the quality vector if quality scores are provided for both sequences.
/// # Arguments
/// - `r1`: A slice of bytes representing the first sequence.
/// - `r1_qual`: An optional slice of bytes representing the quality scores for the first sequence.
/// - `r2`: A slice of bytes representing the second sequence.
/// - `r2_qual`: An optional slice of bytes representing the quality scores for the second sequence.
/// - `overlap`: An `OverlapResult` containing the offset, overlap length, and number of mismatches.
/// # Returns
/// - `EndJoiningResult`: The result of the end joining operation.
///   - If there is no overlap, it returns a concatenated sequence of `r1` and `r2`.
///   - If there is an overlap, it returns a sequence that includes the consensus bases from the overlapping region.
///   - The quality scores are also included if available, with the consensus quality scores calculated based on the overlapping region.
fn join_with_overlap(
    r1: &[u8],
    r1_qual: Option<&[u8]>,
    r2: &[u8],
    r2_qual: Option<&[u8]>,
    overlap: OverlapResult,
) -> EndJoiningResult {
    let offset = overlap.offset;
    let overlap_len = overlap.overlap_len;

    // Calculate prefix regions
    let (prefix_seq, prefix_qual): (Vec<u8>, Option<Vec<u8>>) = if offset < 0 {
        let plen = (-offset) as usize;
        (r2[..plen].to_vec(), r2_qual.map(|q| q[..plen].to_vec()))
    } else {
        let plen = offset as usize;
        (r1[..plen].to_vec(), r1_qual.map(|q| q[..plen].to_vec()))
    };

    // Calculate starting indices for overlap region
    let r1_overlap_start = if offset < 0 { 0 } else { offset as usize };
    let r2_overlap_start = if offset < 0 { (-offset) as usize } else { 0 };

    // Overlap consensus bases and qualities
    let mut overlap_seq = Vec::with_capacity(overlap_len);
    let mut overlap_qual = if r1_qual.is_some() && r2_qual.is_some() {
        Some(Vec::with_capacity(overlap_len))
    } else {
        None
    };

    for i in 0..overlap_len {
        let r1_idx = r1_overlap_start + i;
        let r2_idx = r2_overlap_start + i;
        let base1 = r1[r1_idx];
        let base2 = r2[r2_idx];

        // Determine consensus base
        let consensus_base = if base1 == base2 {
            base1
        } else {
            match (r1_qual, r2_qual) {
                (Some(q1), Some(q2)) => {
                    let q1_val = q1.get(r1_idx).copied().unwrap_or(0);
                    let q2_val = q2.get(r2_idx).copied().unwrap_or(0);
                    if q1_val >= q2_val { base1 } else { base2 }
                }
                _ => b'N',
            }
        };
        overlap_seq.push(consensus_base);

        // Quality: take max if both present
        if let (Some(q1), Some(q2), Some(overlap_q)) = (r1_qual, r2_qual, &mut overlap_qual) {
            let q1_val = q1.get(r1_idx).copied().unwrap_or(0);
            let q2_val = q2.get(r2_idx).copied().unwrap_or(0);
            overlap_q.push(std::cmp::max(q1_val, q2_val));
        }
    }

    // Calculate suffix regions
    let r1_overlap_end = r1_overlap_start + overlap_len;
    let r2_overlap_end = r2_overlap_start + overlap_len;
    let (suffix_seq, suffix_qual): (Vec<u8>, Option<Vec<u8>>) = if r1_overlap_end < r1.len() {
        (
            r1[r1_overlap_end..].to_vec(),
            r1_qual.map(|q| q[r1_overlap_end..].to_vec()),
        )
    } else if r2_overlap_end < r2.len() {
        (
            r2[r2_overlap_end..].to_vec(),
            r2_qual.map(|q| q[r2_overlap_end..].to_vec()),
        )
    } else {
        (vec![], None)
    };

    // Final assembly
    let mut seq = prefix_seq;
    seq.extend_from_slice(&overlap_seq);
    seq.extend_from_slice(&suffix_seq);

    let quality = match (prefix_qual, overlap_qual, suffix_qual) {
        (Some(mut pre), Some(mut ovl), Some(mut suf)) => {
            pre.append(&mut ovl);
            pre.append(&mut suf);
            Some(pre)
        }
        _ => None,
    };

    EndJoiningResult { seq, quality }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_find_best_overlap() {
        let r1 = b"GGGGGGGGGGAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA";
        let r2 = b"AAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATTTTTTTTTT";
        let best_overlap = find_best_overlap(r1, r2, MIN_OVERLAP, 0.0);
        assert_eq!(best_overlap.offset, 10);
        assert_eq!(best_overlap.overlap_len, 100);

        let r1 = b"ACGTACGT";
        let r2 = b"TACGTCG";
        let best_overlap = find_best_overlap(r1, r2, 2, ERROR_RATE_FOR_ENDJOINING);
        assert_eq!(best_overlap.offset, 3);
        assert_eq!(best_overlap.overlap_len, 5);
        assert_eq!(best_overlap.mismatches, 0);

        let r1 = b"GGGGGGGTT";
        let r2 = b"AAAGGGGGGG";
        let best_overlap = find_best_overlap(r1, r2, 2, ERROR_RATE_FOR_ENDJOINING);
        assert_eq!(best_overlap.offset, -3);
        assert_eq!(best_overlap.overlap_len, 7);

        let r1 = b"GGGGGGGTTTTTTTTTTTTTTT";
        let r2 = b"AAAAAAAAAAAAAAGGGGGGG";
        let best_overlap = find_best_overlap(r1, r2, 2, ERROR_RATE_FOR_ENDJOINING);
        assert_eq!(best_overlap.offset, r1.len() as isize);
        assert_eq!(best_overlap.overlap_len, 0);
    }

    #[test]
    fn test_join1() {
        let r1 = b"ACGTACGTTACGT";
        let r2 = b"TACGTTACGTCGA";
        let fasta1 = fasta::Record::with_attrs("r1", None, r1);
        let fasta2 = fasta::Record::with_attrs("r2", None, r2);
        let input = EndJoiningInput::Fasta((&fasta1, &fasta2));
        let overlap = find_best_overlap(r1, r2, MIN_OVERLAP, ERROR_RATE_FOR_ENDJOINING);
        assert_eq!(overlap.offset, 3);
        assert_eq!(overlap.overlap_len, 10);
        let overlap_result =
            OverlapResult::from_simple_overlap(r1.len(), r2.len(), overlap.overlap_len);
        assert_eq!(overlap_result.offset, 3);
        assert_eq!(overlap_result.overlap_len, 10);
        let result = end_joining(
            input.clone(),
            &EndJoiningStrategy::Overlap(overlap.overlap_len),
        )
        .unwrap();
        assert_eq!(result.seq, b"ACGTACGTTACGTCGA");
        let result = end_joining(input.clone(), &EndJoiningStrategy::Overlap(0)).unwrap();
        assert_eq!(result.seq, b"ACGTACGTTACGTTACGTTACGTCGA");
    }

    #[test]
    fn test_join2() {
        let r1 = b"GGGGGGGTT";
        let r2 = b"AAAGGGGGGG";
        let fasta1 = fasta::Record::with_attrs("r1", None, r1);
        let fasta2 = fasta::Record::with_attrs("r2", None, r2);
        let input = EndJoiningInput::Fasta((&fasta1, &fasta2));
        let overlap = find_best_overlap(r1, r2, 4, ERROR_RATE_FOR_ENDJOINING);
        assert_eq!(overlap.offset, -3);
        assert_eq!(overlap.overlap_len, 7);
        let result = join_with_overlap(r1, None, r2, None, overlap.clone());
        assert_eq!(result.seq, b"AAAGGGGGGGTT");

        let result = end_joining(input.clone(), &EndJoiningStrategy::UnknownOverlap);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().seq, b"GGGGGGGTTAAAGGGGGGG");
    }

    #[test]
    fn test_join3() {
        let r1 = b"CCCGGGGGGGTTTTTCCC";
        let r2 = b"GGGGGTTTTTC";

        let fasta1 = fasta::Record::with_attrs("r1", None, r1);
        let fasta2 = fasta::Record::with_attrs("r2", None, r2);
        let input = EndJoiningInput::Fasta((&fasta1, &fasta2));
        let overlap = find_best_overlap(r1, r2, 10, ERROR_RATE_FOR_ENDJOINING);
        assert_eq!(overlap.offset, 5);
        assert_eq!(overlap.overlap_len, 11);
        let result = end_joining(input.clone(), &EndJoiningStrategy::UnknownOverlap);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().seq, b"CCCGGGGGGGTTTTTCCC");
    }

    #[test]
    fn test_join4() {
        let r1 =     b"CAATACATCACAACTGTTTAATAGTACTTGGATTAATGGTACTAGGAAAGGTACTGAAGGAAATGTTACAGAAAATATCATACTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGGAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAAAGGAATGATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGGGATGGTGGTGAGAACAAAAACAAGAGCGAGCCCGAGGTCTTCAGACCTGGAGGAGGAGATATGAGGGACA";
        let r2 =        b"TACATCACAACTGTTTAATAGTACTTGGATTAATGGTACTAGGAAAGGTACTGAAGGAAATGTTACAGAAAATATCATACTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGGAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAAAGGAATGATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGGGATGGTGGTGAGAACAAAAACAAGAGCGAGCCCGAGGTCTTCAGACCTGGAGGAGGAGATATGAGGGAC";
        let joined = b"CAATACATCACAACTGTTTAATAGTACTTGGATTAATGGTACTAGGAAAGGTACTGAAGGAAATGTTACAGAAAATATCATACTCCCATGCAGAATAAAACAAATTATAAACATGTGGCAGGAAGTAGGAAAAGCAATGTATGCCCCTCCCATCAAAGGAATGATTAGATGTTCATCAAATATTACAGGGCTGCTATTAACAAGGGATGGTGGTGAGAACAAAAACAAGAGCGAGCCCGAGGTCTTCAGACCTGGAGGAGGAGATATGAGGGACA";
        let fasta1 = fasta::Record::with_attrs("r1", None, r1);
        let fasta2 = fasta::Record::with_attrs("r2", None, r2);
        let input = EndJoiningInput::Fasta((&fasta1, &fasta2));
        let result = end_joining(input.clone(), &EndJoiningStrategy::UnknownOverlap);
        assert!(result.is_ok());
        assert_eq!(result.unwrap().seq, joined);
    }
}
