use std::collections::HashMap;
use std::error::Error;
use std::fmt::Display;
use std::ops::Range;

use bio::io::fastq::Record;
use getset::{Getters, Setters};
use itertools::Itertools;
use rayon::prelude::*;
use serde::{Deserialize, Serialize};
use virust_locator::prelude::*;

use crate::helper::consensus::{
    self, ConsensusInput, ConsensusParams, ConsensusStrategy, consensus,
};
use crate::helper::end_joining::*;
use crate::helper::params::{QcConfig, TrimConfig};
use crate::helper::tcs_helper::*;
use crate::helper::umis::{UMIDistError, UMIInformationBlocks, UMISummary};

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters)]
pub struct TcsConsensus {
    #[getset(get = "pub", set = "pub")]
    umi_information_block: String,
    #[getset(get = "pub", set = "pub")]
    umi_family_size: usize,
    #[getset(get = "pub", set = "pub")]
    r1_consensus: Record,
    #[getset(get = "pub", set = "pub")]
    r2_consensus: Record,
    #[getset(get = "pub", set = "pub")]
    joined_consensus: Option<Record>,
    #[getset(get = "pub", set = "pub")]
    qc: TcsConsensusQcResult,
    #[getset(get = "pub", set = "pub")]
    trimmed: Option<Record>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Default, PartialEq)]
pub enum TcsConsensusQcResult {
    #[default]
    QcNotInitialized,
    NoJoinedConsensus,
    NotRequired,
    Passed,
    NotPassed(QcNotPassedReport),
    LocatorWithErrors(String),
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Getters, Setters)]
pub struct QcNotPassedReport {
    #[getset(get = "pub")]
    qc_reference: String,
    #[getset(get = "pub")]
    qc_coordinates1: Option<Range<u32>>,
    #[getset(get = "pub")]
    qc_coordinates2: Option<Range<u32>>,
    #[getset(get = "pub")]
    qc_indels: bool,
    #[getset(get = "pub")]
    locator_coordinates: Option<Range<u32>>,
    #[getset(get = "pub")]
    locator_indels: bool,
}

impl Display for QcNotPassedReport {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "QC not passed for reference: {}, coordinates1: {:?}, coordinates2: {:?}, indels: {}, locator_coordinates: {:?}, indels: {}",
            self.qc_reference,
            self.qc_coordinates1,
            self.qc_coordinates2,
            self.qc_indels,
            self.locator_coordinates,
            self.locator_indels,
        )
    }
}

impl Display for TcsConsensusQcResult {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TcsConsensusQcResult::QcNotInitialized => write!(f, "QC not initialized"),
            TcsConsensusQcResult::NoJoinedConsensus => write!(f, "No joined consensus available"),
            TcsConsensusQcResult::NotRequired => write!(f, "QC not required"),
            TcsConsensusQcResult::Passed => write!(f, "QC passed"),
            TcsConsensusQcResult::NotPassed(report) => write!(f, "{}", report),
            TcsConsensusQcResult::LocatorWithErrors(errors) => {
                write!(f, "Locator errors: {}", errors)
            }
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters)]
pub struct TcsConsensusBuildingOutput {
    #[getset(get = "pub")]
    tcs_consensus: Vec<TcsConsensus>,
    #[getset(get = "pub")]
    errors: Vec<String>,
    #[getset(get = "pub")]
    umi_summary: UMISummary,
}

impl TcsConsensus {
    /// Initializes a new `TcsConsensus` instance with empty fields.
    pub fn new() -> Self {
        TcsConsensus {
            umi_information_block: String::new(),
            umi_family_size: 0,
            r1_consensus: Record::new(),
            r2_consensus: Record::new(),
            joined_consensus: None,
            qc: TcsConsensusQcResult::default(),
            trimmed: None,
        }
    }
}

pub fn build_from_filtered_pairs(
    pairs: &Vec<FilteredPair>,
    strategy: consensus::ConsensusStrategy,
    error_cutoff: f32,
) -> Result<TcsConsensusBuildingOutput, UMIDistError> {
    let mut umi_records = HashMap::new();
    let mut umi_information_blocks = Vec::new();
    for pair in pairs {
        let umi_information_block = pair.umi.umi_information_block.clone();
        umi_information_blocks.push(umi_information_block.clone());
        umi_records
            .entry(umi_information_block)
            .or_insert_with(|| Vec::new())
            .push((&pair.r1, &pair.r2));
    }

    let umis = UMIInformationBlocks {
        umi_information_blocks,
    };

    let (umi_families, umi_summary) = umis.find_umi_family_by_error_cutoff(error_cutoff)?;

    let tcs_consensus_results: Vec<Result<TcsConsensus, Box<dyn Error + Send + Sync>>> =
        umi_families
            .families
            .par_iter()
            .map(|umi_family| {
                let umi_information_block = umi_family.umi_information_block.clone();
                let filtered_pairs = umi_records.get(&umi_information_block).ok_or_else(|| {
                    TcsError::UnexpectedError(format!(
                        "No filtered pairs found for UMI information block: {}",
                        umi_information_block
                    ))
                })?;

                let r1_vec = filtered_pairs
                    .iter()
                    .map(|(r1, _)| (*r1).clone())
                    .collect::<Vec<_>>();
                let r2_vec = filtered_pairs
                    .iter()
                    .map(|(_, r2)| (*r2).clone())
                    .collect::<Vec<_>>();

                let r1_consensus =
                    consensus::consensus(strategy, consensus::ConsensusInput::Fastq(&r1_vec))?;
                let r2_consensus =
                    consensus::consensus(strategy, consensus::ConsensusInput::Fastq(&r2_vec))?;

                let r1_consensus_record = Record::with_attrs(
                    &format!("{}_{}_r1", umi_information_block, umi_family.frequency),
                    None,
                    &r1_consensus.seq,
                    &r1_consensus.qual.unwrap(),
                );
                let r2_consensus_record = Record::with_attrs(
                    &format!("{}_{}_r2", umi_information_block, umi_family.frequency),
                    None,
                    &r2_consensus.seq,
                    &r2_consensus.qual.unwrap(),
                );

                let mut tcs_consensus = TcsConsensus::new();
                tcs_consensus.set_umi_information_block(umi_information_block);
                tcs_consensus.set_umi_family_size(umi_family.frequency);
                tcs_consensus.r1_consensus = r1_consensus_record;
                tcs_consensus.r2_consensus = r2_consensus_record;
                tcs_consensus.set_umi_family_size(umi_family.frequency);

                Ok(tcs_consensus)
            })
            .collect();

    let mut tcs_consensus = Vec::new();
    let mut errors = Vec::new();
    for result in tcs_consensus_results {
        match result {
            Ok(consensus) => tcs_consensus.push(consensus),
            Err(e) => errors.push(e.to_string()),
        }
    }
    Ok(TcsConsensusBuildingOutput {
        tcs_consensus,
        errors,
        umi_summary,
    })
}

/// Joins the R1 and R2 consensus FASTQ records into a single joined consensus record.
/// This function takes a mutable reference to a vector of `TcsConsensus` records and performs end joining based on the specified strategy.
/// It will mutate the original `TcsConsensus` records by setting the `joined_consensus` field with the joined record.
/// The joining strategy is determined by the `end_joining_option` parameter.
/// - `tcs_consensus`: A mutable reference to a vector of `TcsConsensus` records.
/// - `end_joining_option`: An integer that determines the joining strategy:
///   - `1`: Simple end joining.
///   - `2`: Overlap end joining with a specified overlap length.
///   - `3`: Unknown overlap, which is not recommended due to potential issues in libraries with many off-target reads.
/// - `overlap_len`: The length of the overlap to use for the overlap end joining strategy.
/// In the orginal Ruby version of TCS, we had an option of 3 for Unknown Overlap but use a consensus strategy to determine the overlap.
/// We often have issues with this approach, particularly in libraries with many off-target reads.
/// We decide to have either with a known overlap length or create end-joining for individual TCS consensus records.
/// The param files generated from previous versions will still be compatible with this function. All option 3 will be converted to EndJoiningStrategy::UnknownOverlap.
///
/// Returns a `Result` indicating success or an error message if joining fails.
/// This function uses parallel processing to join the consensus records efficiently.
/// If any errors occur during the joining process, they are collected and returned as a single error message.
/// If the joining is successful, the joined consensus record is set in the `joined_consensus` field of each `TcsConsensus` record.
pub fn join_consensus_fastq_vec(
    tcs_consensus: &mut Vec<TcsConsensus>,
    end_joining_option: u32,
    overlap_len: usize,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let strategy = match end_joining_option {
        1 => EndJoiningStrategy::Simple,
        2 => EndJoiningStrategy::SimpleOverlap(overlap_len),
        3 => EndJoiningStrategy::Overlap(find_consensus_overlap(
            tcs_consensus
                .iter()
                .map(|c| c.r1_consensus.clone())
                .collect(),
            tcs_consensus
                .iter()
                .map(|c| c.r2_consensus.clone())
                .collect(),
        )?),
        _ => EndJoiningStrategy::UnknownOverlap,
    };

    let errors = tcs_consensus
        .par_iter_mut()
        .filter_map(|consensus| {
            let end_joining_input =
                EndJoiningInput::Fastq((&consensus.r1_consensus, &consensus.r2_consensus));
            let joined_consensus = end_joining(end_joining_input, &strategy);
            match joined_consensus {
                Ok(joined) => {
                    let id = format!(
                        "{}_{}_joined",
                        consensus.umi_information_block, consensus.umi_family_size
                    );
                    let joined_record = Record::with_attrs(
                        &id,
                        None,
                        &joined.seq(),
                        &joined.quality().as_ref().unwrap(),
                    );
                    consensus.set_joined_consensus(Some(joined_record));
                    None
                }
                Err(e) => {
                    consensus.joined_consensus = None;
                    Some(format!(
                        "Error joining consensus for UMI {}: {}",
                        consensus.umi_information_block, e
                    ))
                }
            }
        })
        .collect::<Vec<_>>();
    if !errors.is_empty() {
        return Err(errors.join(", ").into());
    }

    Ok(())
}

const QC_ALGORITHM: QcAlgorithm = QcAlgorithm::SemiGlobal;

pub fn qc_and_trim_consensus_fastq_vec(
    tcs_consensus: &mut Vec<TcsConsensus>,
    qc_config: Option<&QcConfig>,
    trim_config: Option<&TrimConfig>,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    if qc_config.is_none() {
        return Ok(());
    }
    let qc_config = qc_config.unwrap();

    let mut joined_tcs_vec = Vec::new();
    for consensus in tcs_consensus.iter() {
        if let Some(joined) = &consensus.joined_consensus {
            joined_tcs_vec.push(joined.seq());
        }
    }

    let unique_joined_tcs_vec = joined_tcs_vec.into_iter().unique().collect::<Vec<_>>();

    let tcs_qc_input = TcsQcInput::with_attrs(
        unique_joined_tcs_vec,
        qc_config.reference.clone(),
        QC_ALGORITHM,
    )
    .ok_or("Failed to create TcsQcInput")?;

    let qc_output = tcs_qc_input.run_locator()?.results_map().to_owned();

    for consensus in tcs_consensus.iter_mut() {
        if let Some(joined) = &consensus.joined_consensus {
            let joined_seq = joined.seq();
            let joined_qual = joined.qual().to_owned();
            match qc_output.get(joined_seq) {
                Some(Some(locator)) => {
                    let qc_result = get_qc_results(qc_config, locator);
                    consensus.set_qc(qc_result.clone());

                    if trim_config.is_none() {
                        consensus.set_trimmed(None);
                    } else if qc_result == TcsConsensusQcResult::Passed {
                        let trimmed = trim_sequence_from_locator(
                            locator,
                            trim_config.as_ref().unwrap().start as usize,
                            trim_config.as_ref().unwrap().end as usize,
                        );

                        match trimmed {
                            Ok(trimmed) => {
                                consensus.set_trimmed(Some(Record::with_attrs(
                                    &format!(
                                        "{}_{}_trimmed",
                                        consensus.umi_information_block, consensus.umi_family_size
                                    ),
                                    None,
                                    &trimmed.0,
                                    &joined_qual[trimmed.1], // Use the range from the locator to get the quality scores
                                )));
                            }
                            Err(e) => {
                                consensus
                                    .set_qc(TcsConsensusQcResult::LocatorWithErrors(e.to_string()));
                                consensus.set_trimmed(None);
                            }
                        }
                    } else {
                        consensus.set_trimmed(None);
                    }
                }
                Some(None) => {
                    consensus.set_qc(TcsConsensusQcResult::LocatorWithErrors(
                        "Locator returned None".to_string(),
                    ));
                    consensus.set_trimmed(None);
                }
                None => {
                    consensus.set_qc(TcsConsensusQcResult::LocatorWithErrors(format!(
                        "No locator found for sequence: {}",
                        String::from_utf8_lossy(joined_seq)
                    )));
                    consensus.set_trimmed(None);
                }
            }
        } else {
            consensus.set_qc(TcsConsensusQcResult::NoJoinedConsensus);
            consensus.set_trimmed(None);
        }
    }

    Ok(())
}

fn get_qc_results(qc_config: &QcConfig, locator: &Locator) -> TcsConsensusQcResult {
    let locator_ref_start = locator.ref_start;
    let locator_ref_end = locator.ref_end;
    let locator_indel = locator.indel;

    if qc_config.start.is_none() && qc_config.end.is_none() {
        return TcsConsensusQcResult::NotRequired;
    } else if qc_config.start.is_none() {
        if qc_config
            .end
            .as_ref()
            .unwrap()
            .contains(&(locator_ref_end as u32))
            && process_indel_logic(qc_config.indel, locator_indel)
        {
            return TcsConsensusQcResult::Passed;
        }
    } else if qc_config.end.is_none() {
        if qc_config
            .start
            .as_ref()
            .unwrap()
            .contains(&(locator_ref_start as u32))
            && process_indel_logic(qc_config.indel, locator_indel)
        {
            return TcsConsensusQcResult::Passed;
        }
    } else {
        if qc_config
            .start
            .as_ref()
            .unwrap()
            .contains(&(locator_ref_start as u32))
            && qc_config
                .end
                .as_ref()
                .unwrap()
                .contains(&(locator_ref_end as u32))
            && process_indel_logic(qc_config.indel, locator_indel)
        {
            return TcsConsensusQcResult::Passed;
        }
    }

    TcsConsensusQcResult::NotPassed(QcNotPassedReport {
        qc_reference: qc_config.reference.clone(),
        qc_coordinates1: qc_config.start.clone(),
        qc_coordinates2: qc_config.end.clone(),
        qc_indels: qc_config.indel,
        locator_coordinates: Some(locator_ref_start as u32..locator_ref_end as u32),
        locator_indels: locator_indel,
    })
}

pub fn count_passed(tcs_consensus_vec: &Vec<TcsConsensus>) -> usize {
    tcs_consensus_vec
        .iter()
        .filter(|consensus| matches!(consensus.qc, TcsConsensusQcResult::Passed))
        .count()
}

fn process_indel_logic(param_indel_bool: bool, locator_indel: bool) -> bool {
    param_indel_bool || !locator_indel
}

fn find_consensus_overlap(
    r1_consensus: Vec<Record>,
    r2_consensus: Vec<Record>,
) -> Result<OverlapResult, Box<dyn Error + Send + Sync>> {
    let consensus_params = ConsensusParams::default();
    let strategy = ConsensusStrategy::Weighted(consensus_params);
    let r1_consensus_input = ConsensusInput::Fastq(&r1_consensus);
    let r2_consensus_input = ConsensusInput::Fastq(&r2_consensus);
    let r1_consensus_of_consensus = consensus(strategy, r1_consensus_input)?;
    let r2_consensus_of_consensus = consensus(strategy, r2_consensus_input)?;

    Ok(find_best_overlap(
        &r1_consensus_of_consensus.seq,
        &r2_consensus_of_consensus.seq,
        MIN_OVERLAP,
        ERROR_RATE_FOR_ENDJOINING,
    ))
}

#[cfg(test)]
mod tests {

    use super::*;

    #[test]
    fn test_get_qc_results() {
        let qc_config1 = QcConfig {
            reference: "HXB2".to_string(),
            start: Some(6585..6686),
            end: Some(7208..7209),
            indel: true,
        };

        let qc_config2 = QcConfig {
            reference: "HXB2".to_string(),
            start: None,
            end: None,
            indel: true,
        };

        let qc_config3 = QcConfig {
            reference: "HXB2".to_string(),
            start: Some(6585..6686),
            end: Some(7208..7209),
            indel: false,
        };

        let qc_config4 = QcConfig {
            reference: "HXB2".to_string(),
            start: Some(6580..6670),
            end: Some(7208..7209),
            indel: true,
        };

        let qc_config5 = QcConfig {
            reference: "HXB2".to_string(),
            start: Some(6580..6670),
            end: None,
            indel: true,
        };

        let qc_config6 = QcConfig {
            reference: "HXB2".to_string(),
            start: None,
            end: Some(7208..7209),
            indel: true,
        };

        let locator = Locator {
            ref_start: 6585,
            ref_end: 7208,
            indel: true,
            percent_identity: 99.0,
            query_aligned_string: String::new(),
            ref_aligned_string: String::new(),
        };

        let result1 = get_qc_results(&qc_config1, &locator);
        let result2 = get_qc_results(&qc_config2, &locator);
        let result3 = get_qc_results(&qc_config3, &locator);
        let result4 = get_qc_results(&qc_config4, &locator);
        let result5 = get_qc_results(&qc_config5, &locator);
        let result6 = get_qc_results(&qc_config6, &locator);
        assert_eq!(result1, TcsConsensusQcResult::Passed);
        assert_eq!(result2, TcsConsensusQcResult::NotRequired);
        assert_eq!(
            result3,
            TcsConsensusQcResult::NotPassed(QcNotPassedReport {
                qc_reference: "HXB2".to_string(),
                qc_coordinates1: Some(6585..6686),
                qc_coordinates2: Some(7208..7209),
                qc_indels: false,
                locator_coordinates: Some(6585..7208),
                locator_indels: true,
            })
        );
        assert_eq!(result4, TcsConsensusQcResult::Passed);
        assert_eq!(result5, TcsConsensusQcResult::Passed);
        assert_eq!(result6, TcsConsensusQcResult::Passed);
    }
}
