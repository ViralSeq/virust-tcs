use std::error::Error as StdError;
use std::fmt::Display;

use bio::alphabets;
use serde::Deserializer;
use serde::de::Error;
use serde::{Deserialize, Serialize};

use crate::utils::umi::UMI;

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct Params {
    #[serde(deserialize_with = "string_or_number_to_f32")]
    pub platform_error_rate: f32,
    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub platform_format: u32,

    #[serde(alias = "Email", alias = "EMAIL")]
    pub email: Option<String>,

    pub primer_pairs: Vec<RegionParams>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct RegionParams {
    pub region: String,
    pub forward: String,
    pub cdna: String,

    #[serde(deserialize_with = "string_or_number_to_f32")]
    pub majority: f32,
    pub end_join: bool,

    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub end_join_option: u32,
    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub overlap: u32,

    #[serde(alias = "TCS_qc", alias = "TCS_QC")]
    pub tcs_qc: bool,
    pub ref_genome: String,

    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub ref_start: u32,
    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub ref_end: u32,
    pub indel: bool,
    pub trim: bool,
    pub trim_ref: String,
    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub trim_ref_start: u32,
    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub trim_ref_end: u32,
}

impl Display for Params {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{\n")?;
        write!(f, "  platform_error_rate: {},\n", self.platform_error_rate)?;
        write!(f, "  platform_format: {},\n", self.platform_format)?;
        if let Some(email) = &self.email {
            write!(f, "  email: {},\n", email)?;
        } else {
            write!(f, "  email: null,\n")?;
        }
        write!(f, "  primer_pairs: [\n")?;
        for region in &self.primer_pairs {
            write!(f, "    {},\n", region)?;
        }
        write!(f, "  ]\n")?;
        write!(f, "}}")
    }
}

impl Display for RegionParams {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "{{\n")?;
        write!(f, "  region: {},\n", self.region)?;
        write!(f, "  forward: {},\n", self.forward)?;
        write!(f, "  cdna: {},\n", self.cdna)?;
        write!(f, "  majority: {},\n", self.majority)?;
        write!(f, "  end_join: {},\n", self.end_join)?;
        write!(f, "  end_join_option: {},\n", self.end_join_option)?;
        write!(f, "  overlap: {},\n", self.overlap)?;
        write!(f, "  tcs_qc: {},\n", self.tcs_qc)?;
        write!(f, "  ref_genome: {},\n", self.ref_genome)?;
        write!(f, "  ref_start: {},\n", self.ref_start)?;
        write!(f, "  ref_end: {},\n", self.ref_end)?;
        write!(f, "  indel: {},\n", self.indel)?;
        write!(f, "  trim: {},\n", self.trim)?;
        write!(f, "  trim_ref: {},\n", self.trim_ref)?;
        write!(f, "  trim_ref_start: {},\n", self.trim_ref_start)?;
        write!(f, "  trim_ref_end: {}\n", self.trim_ref_end)?;
        write!(f, "}}")
    }
}

impl Params {
    /// Reads a JSON string and converts it into a `Params` struct.
    /// This function uses the `serde_json` library to parse the JSON string.
    /// It ignores extra fields in the JSON string that are not defined in the `Params` struct.
    pub fn from_json_sting(json_str: &str) -> Result<Self, Box<dyn StdError>> {
        let params = serde_json::from_str(json_str)?;
        Ok(params)
    }

    //TODO: write details of the function
    pub fn validate(&self) -> Result<(), String> {
        self.primer_pairs.iter().try_for_each(|region| {
            if region.forward.len() != region.cdna.len() {
                return Err(format!(
                    "Forward and cDNA primers must be the same length: {} vs {}",
                    region.forward.len(),
                    region.cdna.len()
                ));
            }
            if region.ref_start > region.ref_end {
                return Err(format!(
                    "ref_start must be less than ref_end: {} vs {}",
                    region.ref_start, region.ref_end
                ));
            }
            Ok(())
        })?;

        Ok(())
    }
}

fn string_or_number_to_u32<'de, D>(deserializer: D) -> Result<u32, D::Error>
where
    D: Deserializer<'de>,
{
    let val: serde_json::Value = Deserialize::deserialize(deserializer)?;
    match val {
        serde_json::Value::Number(num) => num
            .as_u64()
            .map(|n| n as u32)
            .ok_or_else(|| Error::custom("Invalid number")),
        serde_json::Value::String(s) => {
            if s.trim().is_empty() {
                Ok(0)
            } else {
                Ok(s.parse::<u32>().unwrap_or(0))
            }
        }
        _ => Ok(0),
    }
}

fn string_or_number_to_f32<'de, D>(deserializer: D) -> Result<f32, D::Error>
where
    D: Deserializer<'de>,
{
    let val: serde_json::Value = Deserialize::deserialize(deserializer)?;
    match val {
        serde_json::Value::Number(num) => num
            .as_f64()
            .map(|n| n as f32)
            .ok_or_else(|| Error::custom("Invalid number")),
        serde_json::Value::String(s) => {
            if s.trim().is_empty() {
                Ok(0.0)
            } else {
                Ok(s.parse::<f32>().unwrap_or(0.0))
            }
        }
        _ => Ok(0.0),
    }
}

pub fn validate_cdna_primer(seq: &str) -> Result<(), Box<dyn StdError>> {
    validate_nt_words(seq)?;

    UMI::identify(seq)?;

    Ok(())
}

pub fn validate_nt_words(seq: &str) -> Result<(), Box<dyn StdError>> {
    if seq.is_empty() {
        return Err("Empty sequence".into());
    };
    let alphabet = alphabets::dna::iupac_alphabet();
    if !alphabet.is_word(seq.as_bytes()) {
        return Err(format!("Invalid nucleotide sequence: {}", seq).into());
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    static JSON_STR: &str = r#"
        {
            "platform_error_rate": 0.01,
            "platform_format": 300,
            "email": "shuntaiz@email.unc.edu",
            "primer_pairs": [
                {
                    "region": "RT",
                    "forward": "AAANNNNNGGGGGG",
                    "cdna": "CCCNNNNNNNNNNNNGGGGGGG",
                    "majority": 0.5,
                    "end_join": true,
                    "end_join_option": 2,
                    "overlap": 30,
                    "tcs_qc": true,
                    "ref_genome": "HXB2",
                    "ref_start": 100,
                    "ref_end": 200,
                    "indel": true,
                    "trim": false,
                    "trim_ref": "HXB2",
                    "trim_ref_start": 50,
                    "trim_ref_end": 150
                }
            ]
        }
        "#;

    #[test]
    fn test_params_from_json() {
        let params: Params = serde_json::from_str(JSON_STR).unwrap();
        assert_eq!(params.platform_error_rate, 0.01);
        assert_eq!(params.platform_format, 300);
        assert_eq!(params.email, Some("shuntaiz@email.unc.edu".to_string()));
        assert_eq!(params.primer_pairs.len(), 1);
        assert_eq!(params.primer_pairs[0].region, "RT");
        assert_eq!(params.primer_pairs[0].forward, "AAANNNNNGGGGGG");
        assert_eq!(params.primer_pairs[0].cdna, "CCCNNNNNNNNNNNNGGGGGGG");
        assert_eq!(params.primer_pairs[0].majority, 0.5);
        assert_eq!(params.primer_pairs[0].end_join, true);
        assert_eq!(params.primer_pairs[0].end_join_option, 2);
        assert_eq!(params.primer_pairs[0].overlap, 30);
    }

    #[test]
    fn test_read_json_into_params() {
        let json = std::fs::read_to_string("tests/data/test_params.json").unwrap();

        let params: Params = serde_json::from_str(&json).unwrap();
        assert_eq!(params.platform_error_rate, 0.02);
    }

    #[test]
    fn test_validate_nt_words() {
        let valid_seq = "ATCGNBR";
        let invalid_seq = "ATCGX";

        assert!(validate_nt_words(valid_seq).is_ok());
        assert!(validate_nt_words(invalid_seq).is_err());
    }
}
