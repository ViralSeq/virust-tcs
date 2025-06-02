use std::collections::HashMap;
use std::error::Error;

use bio::io::fasta;
use bio::io::fastq;
use thiserror::Error;

// MARK: ConsensusParams
/// Consensus parameters for the consensus function.
/// The `k` parameter controls the steepness of the logistic curve.
/// The `q0` parameter controls the horizontal shift of the curve.
/// These parameters are used to adjust the confidence level of the consensus base.
#[derive(Debug, Clone, Copy)]
pub struct ConsensusParams {
    k: f64,
    q0: f64,
}

// MARK: Default CosnesusParams
/// Implementing the Default trait for ConsensusParams
/// allows for easy instantiation with default values.
/// The default values are `k = 0.2` and `q0 = 30.0`.
impl Default for ConsensusParams {
    fn default() -> Self {
        Self { k: 0.2, q0: 30.0 }
    }
}

impl ConsensusParams {
    /// Creates a new instance of ConsensusParams with specified k and q0 values.
    /// # Arguments
    /// * `k` - The steepness parameter of the logistic curve.
    /// * `q0` - The midpoint of the logistic curve (controls horizontal shift).
    pub fn new(k: f64, q0: f64) -> Self {
        Self { k, q0 }
    }

    pub fn k(&self) -> f64 {
        self.k
    }
    pub fn q0(&self) -> f64 {
        self.q0
    }
}

// MARK: ConsensusStrategy
/// Enum for consensus method/strategy.
/// The `Weighted` variant uses a logistic function to adjust the confidence level based on quality scores.
/// The `Supermajority` variant uses a super-majority cutoff.
/// The `SimpleMajority` variant uses a simple majority rule.
pub enum ConsensusStrategy {
    Weighted(ConsensusParams),
    Supermajority(f64),
    SimpleMajority,
}

// MARK: ConsensusInput
/// Enum for input type (FASTA or FASTQ).
/// The `Fastq` variant contains a slice of FASTQ records (bio::io::fastq::record).
/// The `Fasta` variant contains a slice of FASTA records (bio::io::fasta::record).
pub enum ConsensusInput<'a> {
    Fastq(&'a [fastq::Record]),
    Fasta(&'a [fasta::Record]),
}

// MARK: ConsensusResult
/// Output struct, always contain a consensus sequence (Vec<u8>), may have Phred+33 encoded per-base quality scores (for FASTQ only).
/// The `quality` field is optional and is only present if the input was FASTQ.
/// The `seq` field contains the consensus sequence.
/// The `qual` field contains the Phred-scaled quality scores.
#[derive(Debug, Clone)]
pub struct ConsensusResult {
    pub seq: Vec<u8>,
    pub qual: Option<Vec<u8>>,
}

// MARK: ConsensusError
#[derive(Error, Debug)]
pub enum ConsensusError {
    #[error(
        "At least 2 records are required to compute a consensus. Current number of records: {0}"
    )]
    InvalidRecordsNumber(usize),
    #[error("All sequences must be the same length")]
    InvalidSequenceLength,
    #[error("Missing quality scores for the Weighted strategy")]
    MissingQualityScores,
}

// MARK: Consensus function
/// Unified consensus function that takes a consensus strategy and input type.
/// This function computes the consensus sequence based on the specified strategy.
/// It can handle both FASTA and FASTQ inputs.
/// # Arguments
/// * `strategy` - The consensus strategy to use.
/// * `input` - The input type, either FASTA or FASTQ.
/// # Returns
/// * `Result<ConsensusResult, Box<dyn Error>>` - A Result containing the consensus sequence and its quality scores as a ConsensusResult struct or an error.
/// # Errors
/// * Returns an error if the input records are empty or if the sequences are not of the same length.
/// # Example 1, consensus with FASTQ input:
/// ```
/// use virust_tcs::helper::consensus::{consensus, ConsensusStrategy, ConsensusInput, ConsensusParams};
/// use bio::io::fastq;
/// use bio::io::fasta;
/// use std::error::Error;
/// use std::str::from_utf8;
/// fn main() -> Result<(), Box<dyn Error>> {
///     let records = vec![
///         fastq::Record::with_attrs("SEQ_ID", None, b"ACGT", b"IIII"),
///         fastq::Record::with_attrs("SEQ_ID", None, b"ACGT", b"IIII"),
///     ];
///     let params = ConsensusParams::default();
///     let input = ConsensusInput::Fastq(&records);
///     let strategy = ConsensusStrategy::Weighted(params);
///     let consensus = consensus(strategy, input)?;  
///     println!("Consensus sequence: {}", from_utf8(&consensus.seq)?);
///     if let Some(qual) = consensus.qual {
///         println!("Consensus quality: {}", from_utf8(&qual)?);
///     }
///     Ok(())
/// }
/// ```
/// # Example 2, consensus with FASTA input, supermajority strategy:
/// ```
/// use virust_tcs::helper::consensus::{consensus, ConsensusStrategy, ConsensusInput};
/// use bio::io::fasta;
/// use std::error::Error;
/// use std::str::from_utf8;
/// fn main() -> Result<(), Box<dyn Error>> {
///     let records = vec![
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///     ];
///     let cutoff = 0.55;
///     let input = ConsensusInput::Fasta(&records);    
///     let strategy = ConsensusStrategy::Supermajority(cutoff);
///     let consensus = consensus(strategy, input)?;
///     // you expect to see "ACGT" as the consensus sequence, if you use cutoff >= 0.6, you will see "ACGN".
///     println!("Consensus sequence: {}", from_utf8(&consensus.seq)?);
///     Ok(())
/// }
/// ```
/// # Example 3, consensus with FASTA input, simple majority strategy:
/// ```
/// use virust_tcs::helper::consensus::{consensus, ConsensusStrategy, ConsensusInput};
/// use bio::io::fasta;
/// use std::error::Error;
/// use std::str::from_utf8;
/// fn main() -> Result<(), Box<dyn Error>> {
///     let records = vec![
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///     ];
///     let input = ConsensusInput::Fasta(&records);
///     let strategy = ConsensusStrategy::SimpleMajority;
///     let consensus = consensus(strategy, input)?;
///     // you expect to see "ACGT" as the consensus sequence.
///     println!("Consensus sequence: {}", from_utf8(&consensus.seq)?);
///     Ok(())
/// }
/// ```
/// # Note
/// * The `Weighted` strategy requires quality scores.
/// * The `Supermajority` and `SimpleMajority` strategies do not require quality scores.
/// * The `Supermajority` strategy requires a cutoff value.
/// * The `SimpleMajority` strategy does not require a cutoff value.
/// * The `ConsensusResult` struct contains the consensus sequence and optional quality scores.
/// * The `ConsensusError` enum contains error variants for different consensus computation errors.
/// * The `ConsensusInput` enum allows for different input types (FASTA or FASTQ).
/// * The `ConsensusStrategy` enum allows for different consensus strategies (Weighted, Supermajority, SimpleMajority).
/// * The `ConsensusParams` struct contains parameters for the consensus computation.
pub fn consensus(
    strategy: ConsensusStrategy,
    input: ConsensusInput,
) -> Result<ConsensusResult, Box<dyn Error>> {
    // Extract sequence and (optionally) qualities by input type
    let (seqs, quals_opt): (Vec<Vec<u8>>, Option<Vec<Vec<u8>>>) = match input {
        ConsensusInput::Fastq(records) => (
            records.iter().map(|r| r.seq().to_vec()).collect(),
            Some(records.iter().map(|r| r.qual().to_vec()).collect()),
        ),
        ConsensusInput::Fasta(records) => {
            (records.iter().map(|r| r.seq().to_vec()).collect(), None)
        }
    };

    let n_records = seqs.len();
    if n_records < 2 {
        return Err(ConsensusError::InvalidRecordsNumber(n_records).into());
    }

    let seq_len = seqs[0].len();

    if !seqs.iter().all(|r| r.len() == seq_len) {
        return Err(ConsensusError::InvalidSequenceLength.into());
    }
    let mut consensus = Vec::with_capacity(seq_len);
    let mut consensus_quals = Vec::with_capacity(seq_len);

    for i in 0..seq_len {
        let bases = seqs.iter().map(|r| r[i]).collect::<Vec<u8>>();

        match &strategy {
            ConsensusStrategy::Weighted(params) => {
                // need qualities for this strategy
                if let Some(quals) = &quals_opt {
                    let col_quals = quals.iter().map(|r| r[i]).collect::<Vec<u8>>();
                    match consensus_base_column_with_quality(
                        &bases, &col_quals, params.k, params.q0,
                    ) {
                        Some((base, qual)) => {
                            consensus.push(base);
                            consensus_quals.push(qual);
                        }
                        None => {
                            consensus.push(b'N');
                            consensus_quals.push(b'!');
                        }
                    }
                } else {
                    return Err(ConsensusError::MissingQualityScores.into());
                }
            }
            ConsensusStrategy::Supermajority(cutoff) => {
                // Ensure cutoff is within valid range, but won't throw an error, force it to be between 0.5 and 1.0
                let cutoff = cutoff.max(0.5);
                let cutoff = cutoff.min(1.0);
                let base = consensus_base_supermajority(&bases, cutoff);
                consensus.push(base);
            }
            ConsensusStrategy::SimpleMajority => {
                let base = consensus_base_simply_majority(&bases);
                consensus.push(base);
            }
        }
    }

    Ok(ConsensusResult {
        seq: consensus,
        qual: match strategy {
            ConsensusStrategy::Weighted(_) => Some(consensus_quals),
            _ => None,
        },
    })
}

// MARK: helper functions
/// Computes a logistic-transformed probability from a Phred quality score.
/// There is a graph in /resources that compares the original Phred quality score vs. logistic-transformed probability with differetn k and q0 values.
///
/// # Arguments
/// * `q` - The Phred quality score (e.g., 20, 30, etc.)
/// * `k` - The steepness parameter of the logistic curve
/// * `q0` - The midpoint of the logistic curve (controls horizontal shift)
///
/// # Returns
/// * A float between 0 and 1 representing the adjusted confidence
pub fn logistic_quality_prob(q: f64, k: f64, q0: f64) -> f64 {
    1.0 / (1.0 + (-k * (q - q0)).exp())
}

/// Converts a Phred quality score into a confidence probability
pub fn phred_quality_prob(q: f64) -> f64 {
    1.0 - 10f64.powf(-q / 10.0)
}

/// Computes consensus base for one column using logistic-transformed quality scores.
pub fn consensus_base_column(bases: &[u8], quals: &[u8], k: f64, q0: f64) -> Option<u8> {
    let mut base_weights: HashMap<u8, f64> = HashMap::new();

    for (&base, &qual_char) in bases.iter().zip(quals.iter()) {
        let q = (qual_char - 33) as f64; // Phred+33 encoding
        let weight = logistic_quality_prob(q, k, q0);
        *base_weights.entry(base).or_insert(0.0) += weight;
    }

    // Find max weight
    let max_weight = base_weights
        .values()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    // Collect all bases with max weight
    let top_bases: Vec<u8> = base_weights
        .into_iter()
        .filter(|&(_, w)| (w - max_weight).abs() < 1e-6)
        .map(|(b, _)| b)
        .collect();

    if top_bases.len() == 1 {
        Some(top_bases[0])
    } else {
        Some(b'N')
    }
}

/// Computes consensus base and its Phred-scaled quality score for one column.
/// Returns (consensus_base, Phred_quality_score).
pub fn consensus_base_column_with_quality(
    bases: &[u8],
    quals: &[u8],
    k: f64,
    q0: f64,
) -> Option<(u8, u8)> {
    let mut base_weights: HashMap<u8, f64> = HashMap::new();

    for (&base, &qual_char) in bases.iter().zip(quals.iter()) {
        let q = (qual_char - 33) as f64; // Phred+33 encoding
        let weight = logistic_quality_prob(q, k, q0);
        *base_weights.entry(base).or_insert(0.0) += weight;
    }

    // Find max weight
    let max_weight = base_weights
        .values()
        .copied()
        .fold(f64::NEG_INFINITY, f64::max);

    // Collect all bases with max weight
    let top_bases: Vec<u8> = base_weights
        .clone()
        .into_iter()
        .filter(|&(_, w)| (w - max_weight).abs() < 1e-6)
        .map(|(b, _)| b)
        .collect();

    if top_bases.len() == 1 {
        let total_weight: f64 = base_weights.values().sum();
        let p_error = if total_weight > 0.0 {
            1.0 - (max_weight / total_weight)
        } else {
            1.0
        };
        let p_error = p_error.max(1e-10); // Avoid log(0)

        let q_consensus = -10.0 * p_error.log10();
        let q_consensus = q_consensus.min(60.0); // Cap at 93 to avoid overflow

        let qual_byte = q_consensus.round() as u8 + 33; // Convert back to Phred+33
        Some((top_bases[0], qual_byte))
    } else {
        Some((b'N', b'!')) // Return 'N' with low quality
    }
}

/// Compute consensus base at a position using a super-majority cutoff.
pub fn consensus_base_supermajority(bases: &[u8], cutoff: f64) -> u8 {
    let mut counts: HashMap<u8, usize> = HashMap::new();
    let total = bases.len();

    for &base in bases {
        *counts.entry(base).or_insert(0) += 1;
    }

    for (base, count) in counts {
        if (count as f64) / (total as f64) > cutoff {
            return base;
        }
    }

    b'N' // Return 'N' if no base passes the threshold
}

pub fn consensus_base_simply_majority(bases: &[u8]) -> u8 {
    let mut counts: HashMap<u8, usize> = HashMap::new();

    for &base in bases {
        *counts.entry(base).or_insert(0) += 1;
    }

    let max_count = counts.values().copied().max().unwrap_or(0);

    let top_bases: Vec<u8> = counts
        .into_iter()
        .filter(|&(_, count)| count == max_count)
        .map(|(b, _)| b)
        .collect();

    if top_bases.len() == 1 {
        top_bases[0]
    } else {
        b'N' // Return 'N' if there's a tie
    }
}

// MARK: Tests
#[cfg(test)]
mod tests {
    use super::*;
    use bio::io::fastq;

    #[test]
    fn test_consensus_fastq() {
        let records = vec![
            fastq::Record::with_attrs("SEQ_ID", None, b"ACGT", b"IIII"),
            fastq::Record::with_attrs("SEQ_ID", None, b"ACGT", b"IIII"),
        ];
        let params = ConsensusParams::default();
        let input = ConsensusInput::Fastq(&records);
        let strategy = ConsensusStrategy::Weighted(params);
        let consensus = consensus(strategy, input).unwrap();
        println!("Consensus results: {:?}", consensus);
        assert_eq!(consensus.seq, b"ACGT");
    }

    #[test]
    fn test_consensus_fastq_panic() {
        let records = vec![fastq::Record::with_attrs("SEQ_ID", None, b"ACGT", b"IIII")];
        let params = ConsensusParams::default();
        let input = ConsensusInput::Fastq(&records);
        let strategy = ConsensusStrategy::Weighted(params);
        let result = consensus(strategy, input);
        assert!(result.is_err());
        if let Err(e) = result {
            assert_eq!(
                e.to_string(),
                "At least 2 records are required to compute a consensus. Current number of records: 1"
            );
        }
    }

    #[test]
    fn test_consensus_fastq_2() {
        let records = vec![
            fastq::Record::with_attrs("SEQ_ID", None, b"A", b"I"),
            fastq::Record::with_attrs("SEQ_ID", None, b"A", b"G"),
            fastq::Record::with_attrs("SEQ_ID", None, b"G", b"A"),
            fastq::Record::with_attrs("SEQ_ID", None, b"G", b":"),
            fastq::Record::with_attrs("SEQ_ID", None, b"G", b"="),
        ];
        let params = ConsensusParams::default();
        let input = ConsensusInput::Fastq(&records);
        let strategy = ConsensusStrategy::Weighted(params);
        let consensus = consensus(strategy, input).unwrap();
        println!("Consensus results: {:?}", consensus);
        assert_eq!(consensus.seq, b"A");
    }

    #[test]
    fn test_consensus_fasta() {
        let records = vec![
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
        ];
        let cutoff = 0.55;
        let input = ConsensusInput::Fasta(&records);
        let strategy = ConsensusStrategy::Supermajority(cutoff);
        let consensus = consensus(strategy, input).unwrap();
        assert_eq!(consensus.seq, b"ACGT");
    }

    #[test]
    fn test_consensus_fasta_n() {
        let records = vec![
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
            fasta::Record::with_attrs("SEQ_ID", None, b"ACGC"),
        ];
        let cutoff = 0.55;
        let input = ConsensusInput::Fasta(&records);
        let strategy = ConsensusStrategy::Supermajority(cutoff);
        let consensus = consensus(strategy, input).unwrap();
        assert_eq!(consensus.seq, b"ACGN");
    }
}
