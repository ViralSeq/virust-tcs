use core::ops::Range;
use rand::prelude::*;
use rand_chacha::ChaCha8Rng;
use regex::Regex;
use serde::Deserialize;
use serde::Serialize;
use std::error;

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub enum UMIType {
    /// Unique Molecular Identifier, regular UMI, for example, "N{11}"
    UMI,
    /// Patterned UMI, for example, "N{3}RYN{3}RYN{3}RYN{3}"
    UMIWithPattern,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct UMI {
    pub umi_type: UMIType,
    pub umi_block: String,
    pub information_index: Vec<u32>,
    pub umi_information_block: String,
}

impl UMI {
    /// Identify UMI from the prototype of the sequence (including Ns as UMI)
    /// # Arguments
    /// * `seq` - The sequence string to identify UMI from.
    /// # Returns
    /// * `Result<(Self, Range<usize>), Box<dyn error::Error>>` - A tuple containing the identified UMI and its range in the sequence.
    /// # Errors
    /// * Returns an error if no UMI is found or if more than one regular UMI is found.
    /// # Example
    /// ```
    /// use virust_tcs::utils::umi::*;
    /// let seq = "AAAAAANNNNNNNNNNGGGGG";
    /// let (umi, range) = UMI::identify(seq).unwrap();
    /// assert_eq!(umi.umi_type, UMIType::UMI);
    /// assert_eq!(umi.umi_block, "NNNNNNNNNN");
    /// assert_eq!(umi.information_index, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
    /// assert_eq!(range.start, 6);
    /// assert_eq!(range.end, 16);
    /// ```
    /// ```
    /// use virust_tcs::utils::umi::*;
    /// let seq = "AAAAAANNNRYNNNRYNNNRYNNNGGGGG";
    /// let (umi, range) = UMI::identify(seq).unwrap();
    /// assert_eq!(umi.umi_type, UMIType::UMIWithPattern);
    /// assert_eq!(umi.umi_block, "NNNRYNNNRYNNNRYNNN");
    /// assert_eq!(umi.information_index, vec![0, 1, 2, 5, 6, 7, 10, 11, 12, 15, 16, 17]);
    /// assert_eq!(range.start, 6);
    /// assert_eq!(range.end, 24);
    /// ```
    /// ```
    /// use virust_tcs::utils::umi::*;
    /// let seq = "AAAAAANNNNNNNNNNGGGGGNNNNNNNNN";
    /// let result = UMI::identify(seq);
    /// assert!(result.is_err());
    /// ```
    /// ```
    /// use virust_tcs::utils::umi::*;
    /// let seq = "AAAAAAGGGGG";
    /// let result = UMI::identify(seq);
    /// assert!(result.is_err());
    /// assert!(result.unwrap_err().to_string() == "No UMI found");
    /// ```
    pub fn identify(seq: &str) -> Result<(Self, Range<usize>), Box<dyn error::Error>> {
        //find UMI, and decide if it is a regular UMI or patterned UMI
        //regular UMI is like a string of Ns of at least 8, for example, "N{11}"
        //patterned UMI is like "N{3}RYN{3}RYN{3}RYN{3}", it has 3-4 blocks of Ns, each with 3-4 Ns, with spacer (any `alphabets::dna::iupac_alphabet()` characters other than N, usually 2-3 characters) in between'
        let umi_regex = Regex::new(r"N{8,}")?;
        let umi_with_pattern_regex = Regex::new(r"(N{3,4}[^N]{2,3}){2,5}N{3,4}")?;
        let matching_regular_umi_count = umi_regex.find_iter(seq).count();
        if matching_regular_umi_count == 1 {
            // regular UMI
            let umi_matching = umi_regex.find(seq).unwrap();
            let umi_block = umi_matching.as_str().to_string();
            let information_index = extract_information_index(&umi_block);
            let umi_information_block = umi_block.clone();
            Ok((
                UMI {
                    umi_type: UMIType::UMI,
                    umi_block,
                    information_index,
                    umi_information_block,
                },
                umi_matching.range(),
            ))
        } else if matching_regular_umi_count > 1 {
            Err("More than one Regular UMI found, and they do not fit the patterned UMI format. Please check your input.".to_string().into())
        } else if umi_with_pattern_regex.is_match(seq) {
            // patterned UMI
            let umi_matching = umi_with_pattern_regex.find(seq).unwrap();
            let umi_block = umi_matching.as_str().to_string();
            let information_index = extract_information_index(&umi_block);

            let umi_information_block: String = information_index
                .iter()
                .map(|&i| umi_block.chars().nth(i as usize).unwrap())
                .collect();

            Ok((
                UMI {
                    umi_type: UMIType::UMIWithPattern,
                    umi_block,
                    information_index,
                    umi_information_block,
                },
                umi_matching.range(),
            ))
        } else {
            Err("No UMI found".to_string().into())
        }
    }

    pub fn generate_regular_umi(umi_length: u32, seed: u64) -> UMI {
        let mut rng = ChaCha8Rng::seed_from_u64(seed);
        let bases = ['A', 'T', 'C', 'G'];
        let mut umi = String::new();
        let mut information_index = Vec::new();

        // Generate a single UMI
        for i in 0..umi_length {
            let base = bases
                .choose(&mut rng)
                .expect("Failed to choose a base")
                .to_owned(); // Choose a random base from the array
            umi.push(base);
            information_index.push(i);
        }
        let umi_information_block = umi.clone();
        UMI {
            umi_type: UMIType::UMI,
            umi_block: umi,
            information_index: information_index,
            umi_information_block,
        }
    }
}

/// algorithm to extract information index from UMI, only count position when it is "N"
/// for example, "N{3}RYN{3}RYN{3}RYN{3}" will return [0, 1, 2, 4, 5, 6, 8, 9, 10, 12, 13, 14]
fn extract_information_index(umi: &str) -> Vec<u32> {
    let mut information_index = Vec::new();
    let mut count = 0;
    for c in umi.chars() {
        if c == 'N' {
            information_index.push(count);
            count += 1;
        } else {
            count += 1;
        }
    }
    information_index
}

#[cfg(test)]
mod tests {
    use super::*;

    static REGULAR_UMI: &str = "AAAAAANNNNNNNNNNGGGGG";
    static PATTERNED_UMI: &str = "AAAAAANNNRYNNNRYNNNRYNNNGGGGG";
    static INVALID_UMI: &str = "AAAAAANNNNNNNNNNGGGGGNNNNNNNNN";
    static NO_UMI_FOUND: &str = "AAAAAAGGGGG";

    #[test]
    fn test_identify_umi() {
        let (umi, range) = UMI::identify(REGULAR_UMI).unwrap();
        assert_eq!(umi.umi_type, UMIType::UMI);
        assert_eq!(umi.umi_block, "NNNNNNNNNN");
        assert_eq!(umi.information_index, vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9]);
        assert_eq!(range.start, 6);
        assert_eq!(range.end, 16);

        let (umi, range) = UMI::identify(PATTERNED_UMI).unwrap();
        assert_eq!(umi.umi_type, UMIType::UMIWithPattern);
        assert_eq!(umi.umi_block, "NNNRYNNNRYNNNRYNNN");
        assert_eq!(
            umi.information_index,
            vec![0, 1, 2, 5, 6, 7, 10, 11, 12, 15, 16, 17]
        );
        assert_eq!(range.start, 6);
        assert_eq!(range.end, 24);

        let result = UMI::identify(INVALID_UMI);
        assert!(result.is_err());

        let result = UMI::identify(NO_UMI_FOUND);
        assert!(result.is_err());
        assert!(result.unwrap_err().to_string() == "No UMI found");
    }

    #[test]
    fn test_generate_regular_umi() {
        let umi = UMI::generate_regular_umi(10, 1);
        assert_eq!(umi.umi_type, UMIType::UMI);
        assert_eq!(umi.information_index.len(), 10);
        assert_eq!(umi.umi_block.len(), 10);
        assert_eq!(umi.information_index, (0..10).collect::<Vec<u32>>());
    }
}
