use std::collections::HashMap;
use std::error::Error;

use bio::io::fastq::Record;
use getset::{Getters, Setters};
use rayon::prelude::*;
use serde::{Deserialize, Serialize};

use crate::helper::consensus;
use crate::helper::tcs_helper::FilteredPair;
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
