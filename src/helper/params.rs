use std::collections::HashMap;
use std::error::Error as StdError;
use std::fmt::Display;
use std::ops::Range;

use bio::alphabets;
use once_cell::sync::Lazy;
use regex::Regex;
use serde::de::Error;
use serde::{Deserialize, Deserializer, Serialize};
use thiserror::Error;

use crate::helper::json::FromJsonString;
use crate::helper::umi::UMI;

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
pub struct ValidatedParams {
    pub primer_pairs: Vec<ValidatedRegionParams>,
}

impl ValidatedParams {
    pub fn get_region_params(&self, region: &str) -> Option<&ValidatedRegionParams> {
        self.primer_pairs.iter().find(|p| p.region == region)
    }
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
    #[serde(default, deserialize_with = "string_or_number_to_u32")]
    pub overlap: u32,

    #[serde(alias = "TCS_qc", alias = "TCS_QC")]
    pub tcs_qc: bool,
    pub ref_genome: String,

    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub ref_start: u32,
    #[serde(default, deserialize_with = "string_or_number_to_option_u32")]
    pub ref_start_lower: Option<u32>,
    #[serde(deserialize_with = "string_or_number_to_u32")]
    pub ref_end: u32,
    #[serde(default, deserialize_with = "string_or_number_to_option_u32")]
    pub ref_end_lower: Option<u32>,
    pub indel: bool,
    pub trim: bool,
    pub trim_ref: Option<String>,
    #[serde(default, deserialize_with = "string_or_number_to_option_u32")]
    pub trim_ref_start: Option<u32>,
    #[serde(default, deserialize_with = "string_or_number_to_option_u32")]
    pub trim_ref_end: Option<u32>,
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

pub static PRESETS: Lazy<HashMap<&'static str, &'static str>> = Lazy::new(|| {
    let mut m = HashMap::new();
    m.insert("v1", include_str!("../../resources/dr_presets/v1.json"));
    m.insert("v2", include_str!("../../resources/dr_presets/v2.json"));
    m.insert("v3", include_str!("../../resources/dr_presets/v3.json"));
    m.insert("v4", include_str!("../../resources/dr_presets/v4.json"));
    m
});

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ValidatedRegionParams {
    pub platform_error_rate: f32,
    pub platform_format: u32,
    pub region: String,
    pub forward_matching: ForwardMatching,
    pub cdna_matching: CDNAMatching,
    pub majority: f32,
    pub end_join: bool,
    pub end_join_option: u32,
    pub overlap: u32,
    pub tcs_qc: bool,
    pub qc_config: Option<QcConfig>,
    pub trim: bool,
    pub trim_config: Option<TrimConfig>,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct QcConfig {
    pub reference: String,
    pub start: Option<Range<u32>>,
    pub end: Option<Range<u32>>,
    pub indel: bool,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct TrimConfig {
    pub reference: String,
    pub start: u32,
    pub end: u32,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct ForwardMatching {
    pub forward: String,
    pub leading_n_number: u32,
    pub bio_forward: String,
}

#[derive(Debug, Clone, Deserialize, Serialize)]
pub struct CDNAMatching {
    pub cdna: String,
    pub umi: UMI,
    pub bio_cdna: String,
}

#[derive(Error, Debug)]
pub enum ParamsValidationError {
    #[error("Platform Error rate out of supported range (0..0.1) {0}")]
    InvalidPlatformErrorRate(f32),
    #[error("Biological primer sequence too short, must be at least 6 characters long")]
    ShortBiologicalPrimer,
    #[error("Empty sequence provided")]
    EmptySequence,
    #[error("Invalid nucleotide word, must use IUPAC alphabets: {0}")]
    InValidNucleotideWord(String),
    #[error("Invalid End Join Option, must be between 1 and 4: {0}")]
    InvalidEndJoinOption(u32),
    #[error("Invalid reference genome cooridinates, start {0}, end {1} not valid")]
    InvalidReferenceGenomeCoordinates(u32, u32),
    #[error(
        "If TCS QC is enabled, reference coordinates must be provided for at least one of the start or end. Cannot be both None(0)"
    )]
    TCSQCReferenceCoordinatesNotProvided,
    #[error(
        "Trimming coodrinates on the reference outside of the boundaries of the qc reference coordinates"
    )]
    TrimmingCoordinatesOutsideQCReference,
    #[error("Trimming reference coordinates must be provided, cannot be None")]
    TCSTrimReferenceCoordinatesNotProvided,
    #[error("Request DR params version {0} not supported, supported versions are {1}")]
    UnsupportedDRParamsVersion(String, String),
    #[error("Failed to parse JSON: {0}")]
    JsonParseError(String),
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
        write!(f, "  trim_ref: {:?},\n", self.trim_ref)?;
        write!(f, "  trim_ref_start: {:?},\n", self.trim_ref_start)?;
        write!(f, "  trim_ref_end: {:?}\n", self.trim_ref_end)?;
        write!(f, "}}")
    }
}

impl Params {
    /// Empty constructor for `Params` struct.
    /// This function initializes a `Params` struct with empty fields.
    /// # Returns
    /// * `Params` - A new instance of `Params` with default values.
    pub fn new() -> Self {
        Params {
            platform_error_rate: 0.0,
            platform_format: 0,
            email: None,
            primer_pairs: Vec::new(),
        }
    }

    /// Validate the parameters in the `Params` struct.
    /// This function checks the validity of the platform error rate, primer sequences,
    /// end join options, reference genome coordinates, and other fields.
    /// If any validation fails, it returns a `ParamsValidationError`.
    /// # Returns
    /// * `Result<ValidatedParams, Box<dyn StdError>>` - A result containing validated parameters or an error.
    /// # Errors
    /// * Returns an error if the platform error rate is out of range, primer sequences are invalid,
    ///   end join options are invalid, or reference genome coordinates are invalid.
    pub fn validate(&self) -> Result<ValidatedParams, Box<dyn StdError>> {
        let platform_error_rate = if self.platform_error_rate > 0.1
            || self.platform_error_rate < 0.0
        {
            return Err(
                ParamsValidationError::InvalidPlatformErrorRate(self.platform_error_rate).into(),
            );
        } else {
            self.platform_error_rate
        };
        let platform_format = self.platform_format;
        let mut validated_primer_pairs = Vec::new();

        for primer_pairs in self.primer_pairs.iter() {
            let forward_matching = validate_forward_primer(&primer_pairs.forward)?;
            let cdna_matching = validate_cdna_primer(&primer_pairs.cdna)?;
            if (1..=4).contains(&primer_pairs.end_join_option) == false {
                return Err(ParamsValidationError::InvalidEndJoinOption(
                    primer_pairs.end_join_option as u32,
                )
                .into());
            }

            let mut ref_genome = String::new();
            let mut ref_start = None;
            let mut ref_end = None;
            let mut trim_ref = String::new();
            let mut trim_ref_start = None;
            let mut trim_ref_end = None;

            if primer_pairs.tcs_qc {
                ref_genome = if ["HXB2", "SIVmm239"].contains(&primer_pairs.ref_genome.as_str()) {
                    primer_pairs.ref_genome.clone()
                } else {
                    "HXB2".to_string()
                };
                ref_start =
                    process_qc_ref_number(primer_pairs.ref_start, primer_pairs.ref_start_lower);
                ref_end = process_qc_ref_number(primer_pairs.ref_end, primer_pairs.ref_end_lower);

                match (ref_start.as_ref(), ref_end.as_ref()) {
                    (Some(start), Some(end)) if start.end >= end.start => {
                        return Err(ParamsValidationError::InvalidReferenceGenomeCoordinates(
                            start.end, end.start,
                        )
                        .into());
                    }
                    (None, None) => {
                        return Err(
                            ParamsValidationError::TCSQCReferenceCoordinatesNotProvided.into()
                        );
                    }
                    _ => {}
                }
            }

            if primer_pairs.trim {
                trim_ref = if ["HXB2", "SIVmm239"].contains(&primer_pairs.ref_genome.as_str()) {
                    primer_pairs.ref_genome.clone()
                } else {
                    "HXB2".to_string()
                };

                trim_ref_start = primer_pairs.trim_ref_start;
                trim_ref_end = primer_pairs.trim_ref_end;

                if trim_ref_start.is_none() || trim_ref_end.is_none() {
                    return Err(
                        ParamsValidationError::TCSTrimReferenceCoordinatesNotProvided.into(),
                    );
                }
                if trim_ref_start.as_ref().unwrap() >= trim_ref_end.as_ref().unwrap() {
                    return Err(ParamsValidationError::InvalidReferenceGenomeCoordinates(
                        *trim_ref_start.as_ref().unwrap(),
                        *trim_ref_end.as_ref().unwrap(),
                    )
                    .into());
                }

                if ref_start.is_some()
                    && ref_end.is_some()
                    && ref_start.as_ref().unwrap().end <= *trim_ref_start.as_ref().unwrap()
                    && ref_end.as_ref().unwrap().start >= *trim_ref_start.as_ref().unwrap()
                {
                } else {
                    return Err(ParamsValidationError::TrimmingCoordinatesOutsideQCReference.into());
                }
            }

            validated_primer_pairs.push(ValidatedRegionParams {
                platform_error_rate,
                platform_format,
                region: primer_pairs.region.clone(),
                forward_matching,
                cdna_matching,
                majority: primer_pairs.majority,
                end_join: primer_pairs.end_join,
                end_join_option: primer_pairs.end_join_option,
                overlap: primer_pairs.overlap,
                tcs_qc: primer_pairs.tcs_qc,
                qc_config: if primer_pairs.tcs_qc {
                    Some(QcConfig {
                        reference: ref_genome,
                        start: ref_start,
                        end: ref_end,
                        indel: primer_pairs.indel,
                    })
                } else {
                    None
                },
                trim: primer_pairs.trim,
                trim_config: if primer_pairs.trim {
                    Some(TrimConfig {
                        reference: trim_ref,
                        start: trim_ref_start.unwrap(),
                        end: trim_ref_end.unwrap(),
                    })
                } else {
                    None
                },
            });
        }

        Ok(ValidatedParams {
            primer_pairs: validated_primer_pairs,
        })
    }

    /// Reads a preset name and returns the corresponding `Params` struct.
    /// This function looks up the preset name in the `PRESETS` map and attempts to parse the JSON string
    /// associated with that preset name into a `Params` struct.
    /// # Arguments
    /// * `present_name` - A string slice representing the name of the preset.
    /// # Returns
    /// * `Result<Params, ParamsValidationError>` - A result containing the `Params` struct if successful,
    ///   or a `ParamsValidationError` if the preset name is not found or if there is an error parsing the JSON.
    pub fn from_preset(present_name: &str) -> Result<Self, ParamsValidationError> {
        let mut all_version_names = PRESETS.keys().cloned().collect::<Vec<_>>();
        all_version_names.sort();
        if let Some(json_str) = PRESETS.get(present_name) {
            Params::from_json_string(json_str)
                .map_err(|e| ParamsValidationError::JsonParseError(e.to_string()))
        } else {
            Err(ParamsValidationError::UnsupportedDRParamsVersion(
                present_name.to_string(),
                all_version_names.join(", "),
            ))
        }
    }
}

pub fn dr_presets_names() -> Vec<&'static str> {
    let mut all_version_names = PRESETS.keys().cloned().collect::<Vec<_>>();
    all_version_names.sort();
    all_version_names
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

fn string_or_number_to_option_u32<'de, D>(deserializer: D) -> Result<Option<u32>, D::Error>
where
    D: serde::Deserializer<'de>,
{
    use serde::de::{self, Unexpected, Visitor};
    use std::fmt;

    struct StringOrNumberVisitor;

    impl<'de> Visitor<'de> for StringOrNumberVisitor {
        type Value = Option<u32>;

        fn expecting(&self, formatter: &mut fmt::Formatter) -> fmt::Result {
            formatter.write_str("an integer, a string representing an integer, or nothing")
        }

        fn visit_u64<E>(self, value: u64) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            Ok(Some(value as u32))
        }

        fn visit_str<E>(self, value: &str) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            if value.trim().is_empty() {
                return Ok(None);
            }
            value.parse().map(Some).map_err(|_| {
                de::Error::invalid_value(Unexpected::Str(value), &"a string representing a u32")
            })
        }

        fn visit_none<E>(self) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            Ok(None)
        }

        fn visit_unit<E>(self) -> Result<Self::Value, E>
        where
            E: de::Error,
        {
            Ok(None)
        }
    }

    deserializer.deserialize_any(StringOrNumberVisitor)
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

pub fn validate_cdna_primer(seq: &str) -> Result<CDNAMatching, Box<dyn StdError>> {
    validate_nt_words(seq)?;

    let (umi, umi_range) = UMI::identify(seq)?;

    let bio_range = umi_range.end..seq.len();

    let bio_cdna = &seq
        .get(bio_range)
        .ok_or_else(|| ParamsValidationError::ShortBiologicalPrimer)?;

    if bio_cdna.len() < 6 {
        return Err(ParamsValidationError::ShortBiologicalPrimer.into());
    }

    Ok(CDNAMatching {
        cdna: seq.to_string(),
        umi: umi,
        bio_cdna: bio_cdna.to_string(),
    })
}

pub fn validate_nt_words(seq: &str) -> Result<(), ParamsValidationError> {
    if seq.is_empty() {
        return Err(ParamsValidationError::EmptySequence);
    };
    let alphabet = alphabets::dna::iupac_alphabet();
    if !alphabet.is_word(seq.as_bytes()) {
        return Err(ParamsValidationError::InValidNucleotideWord(
            seq.to_string(),
        ));
    }
    Ok(())
}

pub fn validate_forward_primer(seq: &str) -> Result<ForwardMatching, Box<dyn StdError>> {
    validate_nt_words(seq)?;

    // Match leading N's (at least 3) at the start of the sequence
    let re = Regex::new(r"(N{3,})(.*)")?;
    let (leading_n_number, bio_forward) = if let Some(caps) = re.captures(seq) {
        (
            caps.get(1).map_or(0, |m| m.as_str().len() as u32),
            caps.get(2)
                .map_or_else(|| seq.to_string(), |m| m.as_str().to_string()),
        )
    } else {
        (0, seq.to_string())
    };

    if bio_forward.len() < 6 {
        return Err(ParamsValidationError::ShortBiologicalPrimer.into());
    }
    Ok(ForwardMatching {
        forward: seq.to_string(),
        leading_n_number,
        bio_forward,
    })
}

fn process_qc_ref_number(n1: u32, n2: Option<u32>) -> Option<Range<u32>> {
    if n1 == 0 {
        return None;
    }
    if let Some(n2) = n2 {
        if n1 < n2 { Some(n1..n2) } else { None }
    } else {
        Some(n1..n1 + 1)
    }
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
                    "overlap": "30",
                    "tcs_qc": true,
                    "ref_genome": "HXB2",
                    "ref_start": 100,
                    "ref_start_lower": 120,
                    "ref_end": 200,
                    "indel": true,
                    "trim": true,
                    "trim_ref": "HXB2",
                    "trim_ref_start": 50,
                    "trim_ref_end": 250
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

    #[test]
    fn test_validate_forward_primer() {
        let valid_seq = "CCGGAANNNATCGGAG";
        let valid_seq2 = "ATGGCAGTAGG";
        let invalid_seq = "NNNATCGX";
        let invalid_seq2 = "NNNATCG";

        assert!(validate_forward_primer(valid_seq).is_ok());
        assert!(validate_forward_primer(valid_seq2).is_ok());
        assert!(validate_forward_primer(invalid_seq).is_err());
        if let Err(e) = validate_forward_primer(invalid_seq) {
            assert_eq!(
                e.to_string(),
                "Invalid nucleotide word, must use IUPAC alphabets: NNNATCGX"
            );
        }

        if let Err(e) = validate_forward_primer(invalid_seq2) {
            assert_eq!(
                e.to_string(),
                "Biological primer sequence too short, must be at least 6 characters long"
            );
        }
    }

    #[test]
    fn test_validate_cdna_primer() {
        let valid_seq = "CCGGAANNNNNNNNATCGGAG";
        let valid_seq2 = "ATGGGAACAGTNNNAGNNNAGNNNAGNNNGCAGTAGG";
        let invalid_seq = "NNNATCGX";
        let invalid_seq2 = "NNNATCG";
        let invalid_seq3 = "AAAAAAAAAANNNNNNNNNNCCCC";

        assert!(validate_cdna_primer(valid_seq).is_ok());
        assert_eq!(
            validate_cdna_primer(valid_seq).unwrap().bio_cdna,
            "ATCGGAG".to_string()
        );
        assert_eq!(
            validate_cdna_primer(valid_seq).unwrap().umi.umi_block,
            "NNNNNNNN".to_string()
        );
        assert!(validate_cdna_primer(valid_seq2).is_ok());
        assert!(validate_cdna_primer(invalid_seq).is_err());
        if let Err(e) = validate_cdna_primer(invalid_seq) {
            assert_eq!(
                e.to_string(),
                "Invalid nucleotide word, must use IUPAC alphabets: NNNATCGX"
            );
        }

        if let Err(e) = validate_cdna_primer(invalid_seq2) {
            assert_eq!(e.to_string(), "No UMI found in cDNA primer: NNNATCG");
        }

        if let Err(e) = validate_cdna_primer(invalid_seq3) {
            assert_eq!(
                e.to_string(),
                "Biological primer sequence too short, must be at least 6 characters long"
            );
        }
    }

    #[test]
    fn test_validate_params() {
        let json = std::fs::read_to_string("tests/data/test_params.json").unwrap();

        let params: Params = serde_json::from_str(&json).unwrap();

        let validated_params = params.validate();

        println!("{:?}", validated_params);
        assert!(validated_params.is_ok());
        let validated_params = validated_params.unwrap();
        assert_eq!(
            validated_params.primer_pairs[0]
                .forward_matching
                .bio_forward,
            "TTATGGGATCAAAGCCTAAAGCCATGTGTA"
        );
        assert_eq!(
            validated_params.primer_pairs[0].cdna_matching.umi.umi_block,
            "N".repeat(9)
        );
    }

    #[test]
    fn test_validate_params_invalid() {
        let params: Params = serde_json::from_str(JSON_STR).unwrap();
        let result = params.validate();
        dbg!(&result);
        assert!(result.is_err());
        assert_eq!(
            result.unwrap_err().to_string(),
            "Trimming coodrinates on the reference outside of the boundaries of the qc reference coordinates"
        );
    }

    #[test]
    fn test_preset_params() {
        let preset_name = ["v1", "v2", "v3", "v4"];

        for name in preset_name.iter() {
            let params = Params::from_preset(name);
            assert!(params.is_ok(), "Failed to load preset: {}", name);
            let params = params.unwrap();
            assert!(
                !params.primer_pairs.is_empty(),
                "No primer pairs found in preset: {}",
                name
            );
            let validated_params = params.validate();
            assert!(
                validated_params.is_ok(),
                "Validation failed for preset: {}",
                name
            );
        }
    }

    #[test]
    fn test_preset_params_invalid() {
        let invalid_preset_name = "invalid_preset";
        let params = Params::from_preset(invalid_preset_name);
        dbg!(&params);
        assert!(params.is_err(), "Expected error for invalid preset");
        if let Err(e) = params {
            assert_eq!(
                e.to_string(),
                format!(
                    "Request DR params version {} not supported, supported versions are v1, v2, v3, v4",
                    invalid_preset_name
                )
            );
        }
    }
}
