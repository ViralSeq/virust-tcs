use std::collections::HashMap;
use std::error::Error;

use bio::io::fastq::Record;
use getset::{Getters, Setters};
use itertools::Itertools;
use rayon::{prelude::*, str};
use serde::{Deserialize, Serialize};

use crate::helper::consensus;
use crate::helper::end_joining::EndJoiningInput;
use crate::helper::end_joining::*;
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
    qc: Option<bool>,
    #[getset(get = "pub", set = "pub")]
    trimmed: Option<Record>,
}

impl TcsConsensus {
    pub fn new() -> Self {
        TcsConsensus {
            umi_information_block: String::new(),
            umi_family_size: 0,
            r1_consensus: Record::new(),
            r2_consensus: Record::new(),
            joined_consensus: None,
            qc: None,
            trimmed: None,
        }
    }

    pub fn build_from_filtered_pairs(
        pairs: &Vec<FilteredPair>,
        strategy: consensus::ConsensusStrategy,
        error_cutoff: f32,
    ) -> Result<(Vec<Self>, Vec<String>, UMISummary), UMIDistError> {
        let mut umi_records = HashMap::new();

        for pair in pairs {
            let umi_information_block = pair.umi.umi_information_block.clone();
            umi_records
                .entry(umi_information_block)
                .or_insert_with(|| Vec::new())
                .push((&pair.r1, &pair.r2));
        }

        let umis = UMIInformationBlocks {
            umi_information_blocks: umi_records.keys().cloned().collect::<Vec<_>>(),
        };

        let (umi_families, umi_summary) = umis.find_umi_family_by_error_cutoff(error_cutoff)?;

        let tcs_consensus_results: Vec<Result<TcsConsensus, Box<dyn Error + Send + Sync>>> =
            umi_families
                .families
                .par_iter()
                .map(|umi_family| {
                    let umi_information_block = umi_family.umi_information_block.clone();
                    let filtered_pairs = umi_records.get(&umi_information_block).unwrap();

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
                        &format!("{}-{}-r1", umi_information_block, umi_family.frequency),
                        None,
                        &r1_consensus.seq,
                        &r1_consensus.qual.unwrap(),
                    );
                    let r2_consensus_record = Record::with_attrs(
                        &format!("{}-{}-r2", umi_information_block, umi_family.frequency),
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
        Ok((tcs_consensus, errors, umi_summary))
    }
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
        2 => EndJoiningStrategy::Overlap(overlap_len),
        // In the orginal Ruby version of TCS, we had an option of 3 for Unknown Overlap but use a consensus strategy to determine the overlap.
        // We often had issues with this approach, particularly in libraries with many off-target reads.
        // We decide to have either with a known overlap length or create end-joining for individual TCS consensus records.
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

pub fn qc_consensus_fastq_vec(
    tcs_consensus: &mut Vec<TcsConsensus>,
    reference: String,
    qc_algorithm: u8,
) -> Result<(), Box<dyn Error + Send + Sync>> {
    let mut joined_tcs_vec = Vec::new();
    for consensus in tcs_consensus.iter() {
        if let Some(joined) = &consensus.joined_consensus {
            joined_tcs_vec.push(joined.seq());
        }
    }

    let unique_joined_tcs_vec = joined_tcs_vec.into_iter().unique().collect::<Vec<_>>();

    let _tcs_qc_input = TcsQcInput::with_attrs(unique_joined_tcs_vec, reference, qc_algorithm);

    Ok(())
}
