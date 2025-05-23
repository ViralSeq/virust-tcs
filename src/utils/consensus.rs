use bio::io::fasta;
use bio::io::fastq;
use std::collections::HashMap;
use std::error::Error;

/// Consensus parameters for the consensus function.
/// The `k` parameter controls the steepness of the logistic curve.
/// The `q0` parameter controls the horizontal shift of the curve.
/// These parameters are used to adjust the confidence level of the consensus base.
#[derive(Debug, Clone, Copy)]
pub struct ConsensusParams {
    k: f64,
    q0: f64,
}

/// Implementing the Default trait for ConsensusParams
/// allows for easy instantiation with default values.
/// The default values are `k = 0.2` and `q0 = 30.0`.
impl Default for ConsensusParams {
    fn default() -> Self {
        Self { k: 0.2, q0: 30.0 }
    }
}

/// Computes the consensus sequence from a list of FASTQ records.
/// The consensus is computed using a logistic function to adjust the confidence level
/// based on the quality scores of the bases.
/// The `k` and `q0` parameters control the steepness and horizontal shift of the logistic curve.
/// The function returns the consensus sequence as a String.
/// # Arguments
/// * `records` - A slice of FASTQ records.
/// * `params` - A reference to a ConsensusParams struct containing the parameters for consensus calculation.
/// # Returns
/// * `Result<String, Box<dyn Error>>` - A Result containing the consensus sequence as a String or an error.
/// # Errors
/// * Returns an error if the input records are empty or if the sequences are not of the same length.
/// # Example
/// ```
/// use virust_tcs::utils::consensus::{consensus_fastq, ConsensusParams};
/// use bio::io::fastq;
/// use std::error::Error;
/// fn main() -> Result<(), Box<dyn Error>> {
///     let records = vec![
///         fastq::Record::with_attrs("SEQ_ID", None, b"ACGT", b"IIII"),
///         fastq::Record::with_attrs("SEQ_ID", None, b"ACGT", b"IIII"),
///     ];
///     let params = ConsensusParams::default();
///     let consensus = consensus_fastq(&records, &params)?;
///     println!("Consensus sequence: {}", consensus);
///     Ok(())
/// }
/// ```

pub fn consensus_fastq(
    records: &[fastq::Record],
    params: &ConsensusParams,
) -> Result<String, Box<dyn Error>> {
    let k = params.k;
    let q0 = params.q0;
    if records.len() < 2 {
        return Err("At least 2 records are required to compute a consensus.".into());
    }

    let seq_len = records[0].seq().len();

    if !records.iter().all(|r| r.seq().len() == seq_len) {
        return Err("All sequences must be the same length".into());
    }

    let mut consensus = Vec::with_capacity(seq_len);

    for i in 0..seq_len {
        let bases: Vec<u8> = records.iter().map(|r| r.seq()[i]).collect();
        let quals: Vec<u8> = records.iter().map(|r| r.qual()[i]).collect();

        match consensus_base_column(&bases, &quals, k, q0) {
            Some(base) => consensus.push(base),
            None => consensus.push(b'N'),
        }
    }

    Ok(String::from_utf8(consensus)?)
}

/// Computes the consensus sequence from a list of FASTA records.
/// Since FASTA records do not have quality scores, the consensus is computed
/// using a super-majority cutoff.
/// The function returns the consensus sequence as a String.
/// # Arguments
/// * `records` - A slice of FASTA records.
/// * `cutoff` - A float representing the super-majority cutoff. For example, 0.6 means that a base must appear in greater than 60% of the sequences to be considered the consensus. Minimum 0.5.
/// # Returns
/// * `Result<String, Box<dyn Error>>` - A Result containing the consensus sequence as a String or an error.
/// # Errors
/// * Returns an error if the input records are empty or if the sequences are not of the same length.
/// # Example
/// ```
/// use virust_tcs::utils::consensus::consensus_fasta;
/// use bio::io::fasta;
/// use std::error::Error;
/// fn main() -> Result<(), Box<dyn Error>> {
///     let records = vec![
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGG"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///         fasta::Record::with_attrs("SEQ_ID", None, b"ACGT"),
///     ];
///     let cutoff = 0.6;
///     let consensus = consensus_fasta(&records, cutoff)?;
///     println!("Consensus sequence: {}", consensus);
///    Ok(())
/// }
/// ```
pub fn consensus_fasta(records: &[fasta::Record], cutoff: f64) -> Result<String, Box<dyn Error>> {
    if records.len() < 2 {
        return Err("At least 2 records are required to compute a consensus.".into());
    }

    // Ensure cutoff is within valid range, but won't throw an error, force it to be between 0.5 and 1.0
    let cutoff = cutoff.max(0.5);
    let cutoff = cutoff.min(1.0);

    let seq_len = records[0].seq().len();

    if !records.iter().all(|r| r.seq().len() == seq_len) {
        return Err("All sequences must be the same length".into());
    }

    let mut consensus = Vec::with_capacity(seq_len);

    for i in 0..seq_len {
        let bases: Vec<u8> = records.iter().map(|r| r.seq()[i]).collect();
        let base = consensus_base_supermajority(&bases, cutoff);
        consensus.push(base);
    }

    Ok(String::from_utf8(consensus)?)
}

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
fn consensus_base_column(bases: &[u8], quals: &[u8], k: f64, q0: f64) -> Option<u8> {
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

/// Compute consensus base at a position using a super-majority cutoff.
fn consensus_base_supermajority(bases: &[u8], cutoff: f64) -> u8 {
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
        let consensus = consensus_fastq(&records, &params).unwrap();
        assert_eq!(consensus, "ACGT");
    }

    #[test]
    fn test_consensus_fastq_panic() {
        let records = vec![fastq::Record::with_attrs("SEQ_ID", None, b"ACGT", b"IIII")];
        let params = ConsensusParams::default();
        let result = consensus_fastq(&records, &params);
        assert!(result.is_err());
        if let Err(e) = result {
            assert_eq!(
                e.to_string(),
                "At least 2 records are required to compute a consensus."
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
        let consensus = consensus_fastq(&records, &params).unwrap();
        assert_eq!(consensus, "A");
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
        let consensus = consensus_fasta(&records, cutoff).unwrap();
        assert_eq!(consensus, "ACGT");
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
        let consensus = consensus_fasta(&records, cutoff).unwrap();
        assert_eq!(consensus, "ACGN");
    }
}
