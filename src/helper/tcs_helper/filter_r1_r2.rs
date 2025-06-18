use std::collections::HashMap;
use std::error::Error;
use std::str::from_utf8;

use bio::bio_types::sequence::SequenceRead;
use bio::io::fastq::Record;
use once_cell::sync::Lazy;
use regex::Regex;
use serde::{Deserialize, Serialize};

use crate::helper::params::{CDNAMatching, ForwardMatching, ValidatedParams};
use crate::helper::tcs_helper::*;
use crate::helper::umi::UMI;

// MARK: FilteredPair
#[derive(Debug, PartialEq, Clone)]
pub struct FilteredPair {
    pub region: String,
    pub umi: UMI,
    pub r1: Record, // r1 read in Record format after trimming and filtering
    pub r2: Record, // r2 read in Record format after trimming and filtering
}

// MARK: PairedRecordFilterResult
#[derive(Debug, PartialEq, Clone)]
pub enum PairedRecordFilterResult {
    Valid(FilteredPair),
    Invalid(FilterPairInvalidReason),
}

// MARK: FilterPairInvalidReason
#[derive(Debug, PartialEq, Clone, Eq, Hash, Deserialize, Serialize)]
pub enum FilterPairInvalidReason {
    InvalidRecords(String),
    GeneralFilterFailed(String),
    R1MatchR2Mismatch(String),
    R2MatchR1Mismatch(String),
    R1R2MatchDifferentRegions(String),
    NoMatch(String),
}

// use std::fmt::{self, Display};
// impl Display for FilterPairInvalidReason {
//     fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
//         match self {
//             FilterPairInvalidReason::InvalidRecords(msg) => write!(f, "Invalid records: {}", msg),
//             FilterPairInvalidReason::GeneralFilterFailed(msg) => {
//                 write!(f, "General filter failed: {}", msg)
//             }
//             FilterPairInvalidReason::R1MatchR2Mismatch(region) => {
//                 write!(f, "R1 matches but R2 does not in region(s): {}", region)
//             }
//             FilterPairInvalidReason::R2MatchR1Mismatch(region) => {
//                 write!(f, "R2 matches but R1 does not in region(s): {}", region)
//             }
//             FilterPairInvalidReason::R1R2MatchDifferentRegions(msg) => {
//                 write!(f, "R1 and R2 match different regions: {}", msg)
//             }
//             FilterPairInvalidReason::NoMatch => write!(f, "No match found"),
//         }
//     }
// }

// MARK: filter_r1_r2_pairs

pub fn filter_r1_r2_pairs(
    r1_record: &Record,
    r2_record: &Record,
    params: &ValidatedParams,
) -> Result<PairedRecordFilterResult, Box<dyn Error + Send + Sync>> {
    // need to benchmark this. not sure how much it slows down the process.
    // Also currently it does not raise an error if the records are not valid, just returns an enum variant Invalid() with the error message.
    // It can put populated in later pipelines for debugging purposes.
    match validate_paired_fastq_record(r1_record, r2_record) {
        Ok(_) => {}
        Err(e) => {
            return Ok(PairedRecordFilterResult::Invalid(
                FilterPairInvalidReason::InvalidRecords(e.to_string()),
            ));
        }
    }

    // check if the platform config is valid, and trim the reads to the platform format
    // raise an error if the platform format is longer than the read length.
    let platform_format = params.primer_pairs[0].platform_format as usize;

    let (r1_trunc, r2_trunc) = if (r1_record.len() >= platform_format)
        && (r2_record.len() >= platform_format)
    {
        (
            &r1_record.get_range(0..(platform_format - 1))?,
            &r2_record.get_range(0..(platform_format - 1))?,
        )
    } else {
        return Err(
            TcsError::InvalidReadLength(platform_format, r1_record.len(), r2_record.len()).into(),
        );
    };

    // General quality filter, check for N content or long homopolymers
    match general_filter(r1_trunc, r2_trunc) {
        GeneralFilterResult::Valid => {}
        GeneralFilterResult::Invalid(msg) => {
            return Ok(PairedRecordFilterResult::Invalid(
                FilterPairInvalidReason::GeneralFilterFailed(msg),
            ));
        }
    }

    // To avoid early return of Invalid(), we will collect regions that do not match
    // and return them at the end if no valid pair is found.
    // This is particularly important when some regions share the same forward or cDNA primer,
    // Early return would prevent checking all regions.
    let mut region_no_matches: HashMap<String, FilterPairInvalidReason> = HashMap::new();

    for region_params in &params.primer_pairs {
        let region = &region_params.region;
        let forward_matching = &region_params.forward_matching;
        let cdna_matching = &region_params.cdna_matching;

        // Check if R1 matches the forward matching config
        let r1_match = match r1_matching(r1_trunc, forward_matching) {
            Ok(Some(record)) => Some(record),
            Ok(None) => None, // No match, continue to next region
            Err(e) => return Err(e),
        };

        let r2_match = match r2_matching(r2_trunc, cdna_matching) {
            Ok((Some(umi), Some(record))) => Some((umi, record)),
            Ok((None, None)) => None,
            Err(e) => return Err(e),
            _ => None,
        };

        if r1_match.is_some() && r2_match.is_some() {
            let (umi, r2_record) = r2_match.unwrap();
            let r1_record = r1_match.unwrap();

            // Create the filtered pair
            let filtered_pair = FilteredPair {
                region: region.clone(),
                umi,
                r1: r1_record,
                r2: reverse_complement(&r2_record), // MARK: reverse compl R2
            };

            return Ok(PairedRecordFilterResult::Valid(filtered_pair));
        } else if r1_match.is_some() && r2_match.is_none() {
            // R1 matches but R2 does not

            region_no_matches.insert(
                region.clone(),
                FilterPairInvalidReason::R1MatchR2Mismatch(region.to_string()),
            );
        } else if r1_match.is_none() && r2_match.is_some() {
            // R2 matches but R1 does not
            region_no_matches.insert(
                region.clone(),
                FilterPairInvalidReason::R2MatchR1Mismatch(region.to_string()),
            );
        } else {
            // Neither R1 nor R2 matches
            region_no_matches.insert(
                region.clone(),
                FilterPairInvalidReason::NoMatch("No match".to_string()),
            );
        }
    }

    let no_match_reason = consolidate_no_match(&region_no_matches)?;

    Ok(PairedRecordFilterResult::Invalid(no_match_reason))
}

// MARK: validate_paired_fastq_record

pub fn validate_paired_fastq_record(
    r1_record: &Record,
    r2_record: &Record,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    if r1_record.is_empty() || r2_record.is_empty() {
        return Err(TcsError::EmptyFastqRecord.into());
    }

    let r1_header = r1_record.id().split_whitespace().next().ok_or_else(|| {
        Box::new(TcsError::InvalidR1Header(r1_record.id().to_string()))
            as Box<dyn Error + Send + Sync>
    })?;
    let r2_header = r2_record.id().split_whitespace().next().ok_or_else(|| {
        Box::new(TcsError::InvalidR2Header(r2_record.id().to_string()))
            as Box<dyn Error + Send + Sync>
    })?;

    if r1_header != r2_header {
        return Err(
            TcsError::R1R2HeaderMismatch(r1_header.to_string(), r2_header.to_string()).into(),
        );
    }

    if r1_record.check().is_err() {
        return Err(TcsError::InvalidR1Record(r1_record.id().to_string()).into());
    }
    if r2_record.check().is_err() {
        return Err(TcsError::InvalidR2Record(r2_record.id().to_string()).into());
    }

    Ok(())
}

// MARK: general_filter

static GENERAL_FILTER_REGEX: Lazy<Regex> = Lazy::new(|| {
    Regex::new(r"N|A{11,}|C{11,}|T{11,}|G{11,}").expect("Failed to compile general filter regex")
});

#[derive(Debug, PartialEq)]
pub enum GeneralFilterResult {
    Valid,
    Invalid(String),
}

fn general_filter(r1_record: &Record, r2_record: &Record) -> GeneralFilterResult {
    let r1_seq = from_utf8(r1_record.seq()).unwrap();
    let r2_seq = from_utf8(r2_record.seq()).unwrap();

    let r1_seq = &r1_seq[4..r1_seq.len()]; // the first 4 bases do not contain any information, it is ok to be Ns. 

    let (r1_match, r2_match) = (
        GENERAL_FILTER_REGEX.is_match(r1_seq),
        GENERAL_FILTER_REGEX.is_match(r2_seq),
    );

    if r1_match && !r2_match {
        GeneralFilterResult::Invalid(format!(
            "General filter (N content or long homopolymers indicating quality issues) failed for R1"
        ))
    } else if r2_match && !r1_match {
        GeneralFilterResult::Invalid(format!(
            "General filter (N content or long homopolymers indicating quality issues) failed for R2"
        ))
    } else if r1_match && r2_match {
        GeneralFilterResult::Invalid(format!(
            "General filter (N content or long homopolymers indicating quality issues) failed for both R1 and R2"
        ))
    } else {
        GeneralFilterResult::Valid
    }
}

// MARK: r1_matching and r2_matching
fn r1_matching(
    r1_record: &Record,
    forward_matching: &ForwardMatching,
) -> Result<Option<Record>, Box<dyn Error + Send + Sync>> {
    let r1_seq = from_utf8(r1_record.seq()).ok().unwrap() as &str;
    let bio_forward = &forward_matching.bio_forward;
    let leading_ns = forward_matching.leading_n_number as usize;

    let trim_start_number = leading_ns + bio_forward.len();
    let primer_region_r1 = r1_seq
        .get(leading_ns..trim_start_number)
        .ok_or_else(|| TcsError::InvalidR1Record(r1_record.id().to_string()))?;

    let diff = diff_by_iupac(primer_region_r1, bio_forward).len();

    if diff < 3 {
        return Ok(Some(r1_record.get_range(trim_start_number..r1_seq.len())?));
    } else {
        return Ok(None);
    }
}

fn r2_matching(
    r2_record: &Record,
    cdna_matching: &CDNAMatching,
) -> Result<(Option<UMI>, Option<Record>), Box<dyn Error + Send + Sync>> {
    let r2_seq = from_utf8(r2_record.seq()).ok().unwrap() as &str;
    let bio_cdna = &cdna_matching.bio_cdna;
    let umi_size = cdna_matching.umi.umi_block.len() as usize;

    let r2_umi_block = r2_seq
        .get(0..umi_size)
        .ok_or_else(|| TcsError::InvalidR2Record(r2_record.id().to_string()))?;

    let mut r2_umi_information_block: Vec<char> = Vec::new();

    for &i in cdna_matching.umi.information_index.iter() {
        let my_char = r2_umi_block.chars().nth(i as usize);

        match my_char {
            Some(c) => r2_umi_information_block.push(c),
            None => {
                return Err(TcsError::InvalidR2Record(r2_record.id().to_string()).into());
            }
        }
    }

    let r2_umi_information_block: String = r2_umi_information_block.iter().collect();
    let r2_umi = UMI {
        umi_type: cdna_matching.umi.umi_type.clone(),
        umi_block: r2_umi_block.to_string(),
        information_index: cdna_matching.umi.information_index.clone(),
        umi_information_block: r2_umi_information_block,
    };

    let trim_start_number = umi_size + bio_cdna.len();
    let primer_region_r2 = r2_seq
        .get(umi_size..trim_start_number)
        .ok_or_else(|| TcsError::InvalidR2Record(r2_record.id().to_string()))?;

    // match ambuiguity codes in primer region
    let diff = diff_by_iupac(primer_region_r2, bio_cdna).len();

    if diff < 3 {
        // UMI identification logic can be added here
        Ok((
            Some(r2_umi),
            Some(r2_record.get_range(trim_start_number..r2_seq.len())?),
        ))
    } else {
        Ok((None, None))
    }
}

// MARK: consolidate_no_match
fn consolidate_no_match(
    region_no_matches: &HashMap<String, FilterPairInvalidReason>,
) -> Result<FilterPairInvalidReason, Box<dyn Error + Send + Sync>> {
    if region_no_matches.is_empty() {
        return Err(TcsError::UnexpectedError(
            "Region Not found in R1 R2 filter results".to_string(),
        )
        .into());
    }

    if region_no_matches.len() == 1 {
        return Ok(region_no_matches.values().next().unwrap().clone());
    }

    if region_no_matches
        .values()
        .all(|v| matches!(v, FilterPairInvalidReason::NoMatch(_)))
    {
        return Ok(FilterPairInvalidReason::NoMatch("No match".to_string()));
    }

    // code for R1R2MatchDifferentRegions
    let mut r1_regions = Vec::new();
    let mut r2_regions = Vec::new();
    for (region, reason) in region_no_matches {
        match reason {
            FilterPairInvalidReason::R1MatchR2Mismatch(_) => r1_regions.push(region.clone()),
            FilterPairInvalidReason::R2MatchR1Mismatch(_) => r2_regions.push(region.clone()),
            _ => {}
        }
    }
    // sort the regions to be consistent when comparing.
    r1_regions.sort();
    r2_regions.sort();
    if !r1_regions.is_empty() && !r2_regions.is_empty() {
        return Ok(FilterPairInvalidReason::R1R2MatchDifferentRegions(format!(
            "R1 regions: {}, R2 regions: {}",
            r1_regions.join(", "),
            r2_regions.join(", ")
        )));
    } else if !r1_regions.is_empty() {
        return Ok(FilterPairInvalidReason::R1MatchR2Mismatch(
            r1_regions.join(", "),
        ));
    } else if !r2_regions.is_empty() {
        return Ok(FilterPairInvalidReason::R2MatchR1Mismatch(
            r2_regions.join(", "),
        ));
    } else {
        return Ok(FilterPairInvalidReason::NoMatch("No match".to_string()));
    }
}

// MARK: Tests
#[cfg(test)]

mod tests {
    use super::*;
    use crate::helper::params::{
        ForwardMatching, ValidatedRegionParams, validate_cdna_primer, validate_forward_primer,
    };
    use crate::helper::umi::{UMI, UMIType};
    use bio::io::fastq::Record;

    #[test]
    fn test_validate_paired_fastq_record_ok() {
        let r1_record = Record::with_attrs("myseq1 1:0:0", None, b"ACGT", b"IIII");
        let r2_record = Record::with_attrs("myseq1 2:0:0", None, b"TGCA", b"IIII");

        assert!(validate_paired_fastq_record(&r1_record, &r2_record).is_ok());
    }

    #[test]
    fn test_validate_paired_fastq_record_should_fail() {
        let r1_record = Record::with_attrs("myseq1 1:0:0", None, b"ACGT", b"IIII");
        let r2_record = Record::with_attrs("myseq3 2:0:0", None, b"TGCA", b"IIII");

        assert!(validate_paired_fastq_record(&r1_record, &r2_record).is_err());
        if let Err(e) = validate_paired_fastq_record(&r1_record, &r2_record) {
            assert_eq!(
                e.to_string(),
                "R1 R2 header mismatch: R1: myseq1, R2: myseq3"
            );
        }

        let r1_record = Record::with_attrs("myseq1 1:0:0", None, b"ACGT", b"III");
        let r2_record = Record::with_attrs("myseq1 2:0:0", None, b"TGCA", b"IIII");

        assert!(validate_paired_fastq_record(&r1_record, &r2_record).is_err());
        if let Err(e) = validate_paired_fastq_record(&r1_record, &r2_record) {
            assert_eq!(e.to_string(), "Invalid R1 record: myseq1 1:0:0");
        }

        let r1_record = Record::with_attrs("", None, b"", b"");
        let r2_record = Record::with_attrs("myseq1 2:0:0", None, b"TGCA", b"IIII");
        assert!(validate_paired_fastq_record(&r1_record, &r2_record).is_err());
        if let Err(e) = validate_paired_fastq_record(&r1_record, &r2_record) {
            assert_eq!(e.to_string(), "Empty fastq record");
        }
    }

    #[test]
    fn test_general_filter_valid() {
        let r1_record = Record::with_attrs("test r1", None, b"AGGCCGA", b"IIIIIII");
        let r2_record = Record::with_attrs("test r2", None, b"TCCAGGA", b"IIIIIII");

        assert_eq!(
            general_filter(&r1_record, &r2_record),
            GeneralFilterResult::Valid
        );
    }

    #[test]
    fn test_general_filter_invalid() {
        let r1_record =
            Record::with_attrs("test r1", None, b"ACAGGGGGGGGGGGGA", b"IIIIIIIIIIIIIIII");
        let r2_record =
            Record::with_attrs("test r2", None, b"GGCTACATCTACTGAC", b"IIIIIIIIIIIIIIII");

        assert_eq!(
            general_filter(&r1_record, &r2_record),
            GeneralFilterResult::Invalid(format!(
                "General filter (N content or long homopolymers indicating quality issues) failed for R1"
            ))
        );

        let r1_record =
            Record::with_attrs("test r1", None, b"ACAGGGGGGGGGGGGA", b"IIIIIIIIIIIIIIII");
        let r2_record =
            Record::with_attrs("test r2", None, b"GGCTACANCTACTGAC", b"IIIIIIIIIIIIIIII");

        assert_eq!(
            general_filter(&r1_record, &r2_record),
            GeneralFilterResult::Invalid(format!(
                "General filter (N content or long homopolymers indicating quality issues) failed for both R1 and R2"
            ))
        );
    }

    #[test]
    fn test_r1_matching() {
        let r1_record = Record::with_attrs(
            "test r1",
            None,
            b"NNNNACGTAGCTAGAAAAAAAAAAAAAAAAAAAAAA",
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        );
        let forward_matching = ForwardMatching {
            forward: "NNNNACGTAGCTAGC".to_string(),
            bio_forward: "ACGTAGCTAG".to_string(),
            leading_n_number: 4,
        };

        let result = r1_matching(&r1_record, &forward_matching);

        assert!(result.is_ok());
        if let Ok(Some(record)) = result {
            assert_eq!(record.seq(), b"AAAAAAAAAAAAAAAAAAAAAA");
        } else {
            panic!("Expected a valid R1 match");
        }

        let r1_record = Record::with_attrs(
            "test r1",
            None,
            b"NNNNATGCAGCTAGAAAAAAAAAAAAAAAAAAAAAA",
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        );

        let result = r1_matching(&r1_record, &forward_matching);

        assert!(result.is_ok());
        if let Ok(Some(record)) = result {
            assert_eq!(record.seq(), b"AAAAAAAAAAAAAAAAAAAAAA");
        } else {
            panic!("Expected a valid R1 match");
        }

        let r1_record = Record::with_attrs(
            "test r1",
            None,
            b"NNNNATGCGCCTAGAAAAAAAAAAAAAAAAAAAAAA",
            b"IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII",
        );

        let result = r1_matching(&r1_record, &forward_matching);

        assert!(result.is_ok());
        assert!(result.unwrap().is_none(), "Expected no match for R1");
    }

    #[test]
    fn test_r2_matching() {
        let cdna_primer =
            "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNNNCAGTCCATTTTGCTYTAYTRABVTTACAATRTGC";
        let cdna_matching = validate_cdna_primer(cdna_primer).unwrap();

        let r2_record = Record::with_attrs(
            "test r2",
            None,
            b"TACTGTTTTACCAGTCCATTTTGCTCTATTGACGTTACAATGTGCTTGTCTCATATTTCCTATTTTTCCTATTGTAACAAATGCTCTCCCTGGTCCCCTCTGGATACGGATACTTTTTCTTGTATTGTTGTTGGGTCTTGTACAATTAATTTCTACAGATGTGTTCAGCTGTACTATTATGGTTTTAGCATTGTCCGTGAAATTGACAGATCTAATTACTACCTCTTCTTCTGCTAGACTGCCATTTAACAGCAGTTGAGTTGATACTACTGGCCTAATTCCATGTGTACATTGTACTGT",
            b"CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGEFCGGGFGGGFFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG9BEFGDGGGGGGGGGGGGGGGGGGFGGGGGFGGGGGGGGFFEFGGGGGGGGGFFFGGGGGGGFGFAAFFCGGGGGGGGGCFFGGGGGGGGGGEDGFGGGFGGGGGDFFFFGGGGCFFGGF8DGGGGFGGGGGFF<DBFFGFEEFFGGGFFFFFCEFEEFFFFFFFFFFEEF9@DECEEFEEEECE?EEFFFECEF4*");
        let result = r2_matching(&r2_record, &cdna_matching);
        assert!(result.is_ok());
        assert_eq!(
            result.unwrap(),
            (
                Some(UMI {
                    umi_type: UMIType::UMI,
                    umi_block: "TACTGTTTTAC".to_string(),
                    information_index: vec![0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10],
                    umi_information_block: "TACTGTTTTAC".to_string(),
                }),
                Some(Record::with_attrs(
                    "test r2",
                    None,
                    b"TTGTCTCATATTTCCTATTTTTCCTATTGTAACAAATGCTCTCCCTGGTCCCCTCTGGATACGGATACTTTTTCTTGTATTGTTGTTGGGTCTTGTACAATTAATTTCTACAGATGTGTTCAGCTGTACTATTATGGTTTTAGCATTGTCCGTGAAATTGACAGATCTAATTACTACCTCTTCTTCTGCTAGACTGCCATTTAACAGCAGTTGAGTTGATACTACTGGCCTAATTCCATGTGTACATTGTACTGT",
                    b"GGGGGEFCGGGFGGGFFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG9BEFGDGGGGGGGGGGGGGGGGGGFGGGGGFGGGGGGGGFFEFGGGGGGGGGFFFGGGGGGGFGFAAFFCGGGGGGGGGCFFGGGGGGGGGGEDGFGGGFGGGGGDFFFFGGGGCFFGGF8DGGGGFGGGGGFF<DBFFGFEEFFGGGFFFFFCEFEEFFFFFFFFFFEEF9@DECEEFEEEECE?EEFFFECEF4*",
                )),
            ),
        );
    }

    #[test]
    fn test_filter_r1_r2_pairs() {
        let r2_record = Record::with_attrs(
            "M01825:522:000000000-C7M6N:1:1101:13543:1027 2:N:0:GCCTTAA",
            None,
            b"TACTGTTTTACCAGTCCATTTTGCTCTATTGACGTTACAATGTGCTTGTCTCATATTTCCTATTTTTCCTATTGTAACAAATGCTCTCCCTGGTCCCCTCTGGATACGGATACTTTTTCTTGTATTGTTGTTGGGTCTTGTACAATTAATTTCTACAGATGTGTTCAGCTGTACTATTATGGTTTTAGCATTGTCCGTGAAATTGACAGATCTAATTACTACCTCTTCTTCTGCTAGACTGCCATTTAACAGCAGTTGAGTTGATACTACTGGCCTAATTCCATGTGTACATTGTACTGT",
            b"CCCCCGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGEFCGGGFGGGFFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGG9BEFGDGGGGGGGGGGGGGGGGGGFGGGGGFGGGGGGGGFFEFGGGGGGGGGFFFGGGGGGGFGFAAFFCGGGGGGGGGCFFGGGGGGGGGGEDGFGGGFGGGGGDFFFFGGGGCFFGGF8DGGGGFGGGGGFF<DBFFGFEEFFGGGFFFFFCEFEEFFFFFFFFFFEEF9@DECEEFEEEECE?EEFFFECEF4*");

        let r1_record = Record::with_attrs(
            "M01825:522:000000000-C7M6N:1:1101:13543:1027 1:N:0:GCCTTAA",
            None,
            b"NGAGTTATGGGATCAAAGCCTAAAGCCATGTGTAAAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGGGGAATGCTACTAATACCAATAGTAGTAATACCAATAGTAGTAGCGGGGAAATGATGATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAACATAAGAGGTAAGGTGCAGAAAGAATATGCATTTTTTTATAAACTTGATATAGTACCAATAGATAATACCAGCTATAGGTTGATAAGTTGTAACATCTCAGTCATTACACAGGCCTGTCC",
            b"#8ACCGGGFGG9FEFGGGGGGGEGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGFFGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGGDGGGGGGGFGGGGGGGGGGGGGGGGFFGGGGGGGGGGGDGGCGGGGGGFGGFGGGGGFGGF=CFFFFFCFFFFEEAAFFEEF;D6EFE8;",
        );

        let cdna_primer =
            "GTGACTGGAGTTCAGACGTGTGCTCTTCCGATCTNNNNNNNNNNNCAGTCCATTTTGCTYTAYTRABVTTACAATRTGC";
        let cdna_matching = validate_cdna_primer(cdna_primer).unwrap();

        let forward_primer =
            "GCCTCCCTCGCGCCATCAGAGATGTGTATAAGAGACAGNNNNTTATGGGATCAAAGCCTAAAGCCATGTGTA";
        let forward_matching = validate_forward_primer(forward_primer).unwrap();

        let region_params = ValidatedRegionParams {
            platform_error_rate: 0.01,
            platform_format: 300,
            region: "test_region".to_string(),
            forward_matching,
            cdna_matching,
            majority: 0.6,
            end_join: false,
            end_join_option: 1,
            overlap: 0,
            tcs_qc: false,
            qc_config: None,
            trim: false,
            trim_config: None,
        };

        let validated_params = ValidatedParams {
            primer_pairs: vec![region_params],
        };

        let result = filter_r1_r2_pairs(&r1_record, &r2_record, &validated_params);
        assert!(result.is_ok());
        let paired_result = result.unwrap();
        dbg!(&paired_result);
        match paired_result {
            PairedRecordFilterResult::Valid(filtered_pair) => {
                assert_eq!(filtered_pair.region, "test_region");
                assert_eq!(filtered_pair.umi.umi_information_block, "TACTGTTTTAC");
                assert_eq!(filtered_pair.r1.seq(), b"AAATTAACCCCACTCTGTGTTAGTTTAAAGTGCACTGATTTGGGGAATGCTACTAATACCAATAGTAGTAATACCAATAGTAGTAGCGGGGAAATGATGATGGAGAAAGGAGAGATAAAAAACTGCTCTTTCAATATCAGCACAAACATAAGAGGTAAGGTGCAGAAAGAATATGCATTTTTTTATAAACTTGATATAGTACCAATAGATAATACCAGCTATAGGTTGATAAGTTGTAACATCTCAGTCATTACACAGGCCTGTC");
                assert_eq!(reverse_complement(&filtered_pair.r2).seq(), b"TTGTCTCATATTTCCTATTTTTCCTATTGTAACAAATGCTCTCCCTGGTCCCCTCTGGATACGGATACTTTTTCTTGTATTGTTGTTGGGTCTTGTACAATTAATTTCTACAGATGTGTTCAGCTGTACTATTATGGTTTTAGCATTGTCCGTGAAATTGACAGATCTAATTACTACCTCTTCTTCTGCTAGACTGCCATTTAACAGCAGTTGAGTTGATACTACTGGCCTAATTCCATGTGTACATTGTACTG");
            }
            PairedRecordFilterResult::Invalid(msg) => panic!("Expected valid pair, got: {:?}", msg),
        }
    }

    #[test]
    fn test_consolidate_no_match() {
        let mut region_no_matches1 = HashMap::new();
        region_no_matches1.insert(
            "RT".to_string(),
            FilterPairInvalidReason::NoMatch("No match".to_string()),
        );
        region_no_matches1.insert(
            "V1V3".to_string(),
            FilterPairInvalidReason::NoMatch("No match".to_string()),
        );

        let mut region_no_matches2 = HashMap::new();
        region_no_matches2.insert(
            "RT".to_string(),
            FilterPairInvalidReason::NoMatch("No match".to_string()),
        );
        region_no_matches2.insert(
            "V1V3".to_string(),
            FilterPairInvalidReason::R1MatchR2Mismatch(
                "R1 Match but R2 not match for V1V3".to_string(),
            ),
        );

        let mut region_no_matches3 = HashMap::new();
        region_no_matches3.insert(
            "RT".to_string(),
            FilterPairInvalidReason::NoMatch("No match".to_string()),
        );
        region_no_matches3.insert(
            "V1V3".to_string(),
            FilterPairInvalidReason::R1MatchR2Mismatch(
                "R1 Match but R2 not match for V1V3".to_string(),
            ),
        );
        region_no_matches3.insert(
            "V1V2".to_string(),
            FilterPairInvalidReason::R1MatchR2Mismatch(
                "R1 Match but R2 not match for V1V2".to_string(),
            ),
        );

        let mut region_no_matches4 = HashMap::new();
        region_no_matches4.insert(
            "RT".to_string(),
            FilterPairInvalidReason::R1MatchR2Mismatch(
                "R1 Match but R2 not match for RT".to_string(),
            ),
        );
        region_no_matches4.insert(
            "V1V3".to_string(),
            FilterPairInvalidReason::R2MatchR1Mismatch(
                "R2 Match but R1 not match for V1V3".to_string(),
            ),
        );

        let result = consolidate_no_match(&region_no_matches1);
        assert!(result.is_ok());
        assert_eq!(
            result.unwrap(),
            FilterPairInvalidReason::NoMatch("No match".to_string())
        );

        let result = consolidate_no_match(&region_no_matches2);
        assert!(result.is_ok());
        assert_eq!(
            result.unwrap(),
            FilterPairInvalidReason::R1MatchR2Mismatch("V1V3".to_string())
        );

        let result = consolidate_no_match(&region_no_matches3);
        assert!(result.is_ok());
        assert_eq!(
            result.unwrap(),
            FilterPairInvalidReason::R1MatchR2Mismatch("V1V2, V1V3".to_string())
        );

        let result = consolidate_no_match(&region_no_matches4);
        assert!(result.is_ok());
        assert_eq!(
            result.unwrap(),
            FilterPairInvalidReason::R1R2MatchDifferentRegions(
                "R1 regions: RT, R2 regions: V1V3".to_string()
            )
        );
    }
}
