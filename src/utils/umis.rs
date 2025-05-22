use crate::utils::umi::UMI;
use itertools::Itertools;
use serde::Deserialize;
use serde::Serialize;
use std::error::Error;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct UMIs {
    pub umis: Vec<UMI>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct UMIFamily {
    pub umi_information_block: String,
    pub frequency: usize,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct UMIFamilies {
    pub families: Vec<UMIFamily>,
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
    /// * `Result<Vec<&str>, Box<dyn Error>>` - A vector of UMI families that exceed the error cutoff.

    pub fn find_umi_family_by_error_cutoff(
        &self,
        error_cutoff: f32,
    ) -> Result<UMIFamilies, Box<dyn Error>> {
        let mut umis: Vec<&str> = self.get_information_blocks();
        let mut families: Vec<UMIFamily> = Vec::new();

        umis.sort();

        let umi_freq = umis.iter().dedup_with_count();

        let mut freq_count: Vec<usize> = Vec::new();
        umi_freq.clone().into_iter().for_each(|(count, _)| {
            freq_count.push(count);
        });

        let max_freq = freq_count.iter().max().unwrap();

        let umi_cut_off = umi_cut_off(*max_freq, Some(error_cutoff));

        umi_freq.into_iter().for_each(|(count, umi)| {
            if count > umi_cut_off as usize {
                let umi = UMIFamily {
                    umi_information_block: umi.to_string(),
                    frequency: count,
                };
                families.push(umi);
            }
        });

        Ok(UMIFamilies { families })
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
