use std::collections::HashMap;

use getset::{Getters, Setters};
use itertools::Itertools;
use serde::Deserialize;
use serde::Serialize;
use thiserror::Error;

use crate::helper::umi::UMI;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct UMIs {
    pub umis: Vec<UMI>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct UMIFamily {
    pub umi_information_block: String,
    pub frequency: usize,
}

impl UMIFamily {
    pub fn to_hash(&self) -> HashMap<String, usize> {
        let mut hash = HashMap::new();
        hash.insert(self.umi_information_block.clone(), self.frequency);
        hash
    }
}

#[derive(Debug, Clone, Deserialize, Serialize, Getters, Setters)]
pub struct UMISummary {
    #[getset(get = "pub")]
    umi_cut_off: usize,
    #[getset(get = "pub")]
    umi_freq: HashMap<String, usize>,
    #[getset(get = "pub")]
    umi_freq_distribution: HashMap<usize, usize>,
}

impl UMISummary {
    pub fn get_passed_umis_hashmap(&self) -> HashMap<String, usize> {
        self.umi_freq()
            .iter()
            .filter(|&(_, &count)| count >= *self.umi_cut_off())
            .map(|(umi, &count)| (umi.clone(), count))
            .collect()
    }
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct UMIFamilies {
    pub families: Vec<UMIFamily>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct UMIInformationBlocks {
    pub umi_information_blocks: Vec<String>,
}

#[derive(Debug, Error)]
pub enum UMIDistError {
    #[error("Total UMIs less than 5")]
    TooFewUMIs,
    #[error("Too few Records")]
    TooFewRecords,
}

impl UMIInformationBlocks {
    pub fn umi_information_blocks(&self) -> Vec<&str> {
        self.umi_information_blocks
            .iter()
            .map(|s| s.as_str())
            .collect()
    }

    pub fn from_umis(umis: &UMIs) -> Self {
        let umi_information_blocks: Vec<String> = umis
            .get_information_blocks()
            .iter()
            .map(|s| s.to_string())
            .collect();
        UMIInformationBlocks {
            umi_information_blocks,
        }
    }

    pub fn find_umi_family_by_error_cutoff(
        &self,
        error_cutoff: f32,
    ) -> Result<(UMIFamilies, UMISummary), UMIDistError> {
        let mut umis: Vec<&str> = self.umi_information_blocks();

        if umis.len() < 5 {
            return Err(UMIDistError::TooFewRecords);
        }

        let mut families: Vec<UMIFamily> = Vec::new();

        umis.sort();

        let umi_freq = umis.iter().dedup_with_count();

        let mut freq_count: Vec<usize> = Vec::new();
        umi_freq.clone().into_iter().for_each(|(count, _)| {
            freq_count.push(count);
        });

        if freq_count.len() < 5 {
            return Err(UMIDistError::TooFewUMIs);
        }

        let mut freq_count_distribution: HashMap<usize, usize> = HashMap::new();
        freq_count.iter().for_each(|&count| {
            *freq_count_distribution.entry(count).or_insert(0) += 1;
        });

        freq_count.sort();
        freq_count.reverse();

        // the latest RUBY version of TCS implementation uses the first 5 frequencies to calculate the cut-off.
        // The Rust version will be consistent with that.
        let max_freq: usize =
            (freq_count[0..5].iter().sum::<usize>() as f64 / 5.0).round() as usize;

        let umi_cut_off = umi_cut_off(max_freq, Some(error_cutoff));

        let mut umi_distribution: HashMap<String, usize> = HashMap::new();
        umi_freq.into_iter().for_each(|(count, umi)| {
            umi_distribution.insert(umi.to_string(), count);
            if count > umi_cut_off as usize {
                let umi = UMIFamily {
                    umi_information_block: umi.to_string(),
                    frequency: count,
                };
                families.push(umi);
            }
        });

        Ok((
            UMIFamilies { families },
            UMISummary {
                umi_cut_off,
                umi_freq: umi_distribution,
                umi_freq_distribution: freq_count_distribution,
            },
        ))
    }
}

impl UMIs {
    /// populate the information blocks of UMIs into a vector of &str.
    pub fn get_information_blocks(&self) -> Vec<&str> {
        let mut information_blocks: Vec<&str> = Vec::new();
        for umi in &self.umis {
            information_blocks.push(&umi.umi_information_block);
        }
        information_blocks
    }

    /// Finds UMI families based on their frequency in the input list.
    ///
    /// # Arguments
    ///
    /// * `error_cutoff` - A float representing the error cutoff for UMI family identification.
    ///
    /// # Returns
    ///
    /// * `Result<UMIFamilies, Box<dyn Error>>` - A vector of UMI families that exceed the error cutoff.

    pub fn find_umi_family_by_error_cutoff(
        &self,
        error_cutoff: f32,
    ) -> Result<(UMIFamilies, UMISummary), UMIDistError> {
        let umi_information_blocks = UMIInformationBlocks::from_umis(self);
        umi_information_blocks.find_umi_family_by_error_cutoff(error_cutoff)
    }
}

/// Calculates the UMI cut-off based on the maximum frequency and error cutoff.
/// The cut-off is determined using a polynomial regression model.
/// # Arguments
/// * `m` - The maximum frequency.
/// * `error_cutoff` - An optional float representing the error cutoff.
/// # Returns
/// * `usize` - The calculated UMI cut-off value.
pub fn umi_cut_off(m: usize, error_cutoff: Option<f32>) -> usize {
    let error_cutoff = (error_cutoff.unwrap_or(0.02) * 1000.0).round() as usize;

    fn poly(coeffs: &[f64], m: usize) -> f64 {
        coeffs
            .iter()
            .enumerate()
            .map(|(i, &c)| c * (m as f64).powi((coeffs.len() - 1 - i) as i32))
            .sum()
    }

    let (coeffs, min_val, upper_bound) = match error_cutoff {
        0..=4 => (
            &[
                -9.59e-27, 3.27e-21, -3.05e-16, 1.2e-11, -2.19e-7, 0.004044, 2.273,
            ][..],
            2,
            None,
        ),
        5..=14 => (
            &[
                1.09e-26, 7.82e-22, -1.93e-16, 1.01e-11, -2.31e-7, 0.00645, 2.872,
            ][..],
            2,
            None,
        ),
        _ => (
            &[
                -1.24e-21, 3.53e-17, -3.90e-13, 2.12e-9, -6.06e-6, 1.80e-2, 3.15,
            ][..],
            2,
            Some(8500),
        ),
    };

    if m <= 10 {
        return 2;
    }

    if let Some(upper) = upper_bound {
        if m > upper {
            return (0.0079 * (m as f64) + 9.4869).round() as usize;
        }
    }

    let n = poly(coeffs, m).round() as usize;
    n.max(min_val)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_umi_cut_off() {
        assert_eq!(umi_cut_off(1000, Some(0.015)), 17);
        assert_eq!(umi_cut_off(1000, Some(0.005)), 9);
        assert_eq!(umi_cut_off(10, None), 2);
        assert_eq!(umi_cut_off(8500, None), 83);
        assert_eq!(umi_cut_off(10000, None), 88);
        assert_eq!(umi_cut_off(10000, Some(0.005)), 53);
        assert_eq!(umi_cut_off(10000, Some(0.001)), 30);
    }
}
