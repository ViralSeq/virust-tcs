use std::collections::HashMap;
use std::fmt::Display;

use chrono::{DateTime, Local};
use getset::{Getters, Setters};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

use crate::helper::params::Params;
use crate::helper::tcs_helper::LOW_ABUNDANCE_THRESHOLD_FOR_RAW_READS;
use crate::helper::tcs_helper::TcsConsensus;
use crate::helper::tcs_helper::filter_r1_r2::FilterPairInvalidReason;
use crate::helper::umis::UMISummary;

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters)]
pub struct TcsReport {
    #[getset(get = "pub", set = "pub")]
    process_start_time: DateTime<Local>,
    #[getset(get = "pub", set = "pub")]
    current_version: String,
    #[getset(get = "pub", set = "pub")]
    input_directory: String,
    #[getset(get = "pub", set = "pub")]
    advanced_settings: AdvancedSettings,

    #[getset(get = "pub", set = "pub")]
    input_params: Params,
    #[getset(get = "pub", set = "pub")]
    total_reads: usize,

    #[getset(get = "pub", set = "pub")]
    failed_match_reasons: Vec<FilterPairInvalidReason>,

    #[getset(get = "pub", set = "pub")]
    region_reports: Vec<RegionReport>,
    #[getset(get = "pub")]
    errors: Vec<String>,
    #[getset(get = "pub")]
    warnings: Vec<TcsReportWarnings>,
    #[getset(get = "pub", set = "pub")]
    process_end_time: DateTime<Local>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters)]
pub struct RegionReport {
    #[getset(get = "pub", set = "pub")]
    region_name: String,
    #[getset(get = "pub", set = "pub")]
    filtered_reads_for_region: usize,
    #[getset(get = "pub", set = "pub")]
    tcs_consensus_results: Option<Vec<TcsConsensus>>,
    #[getset(get = "pub", set = "pub")]
    umi_summary: Option<UMISummary>,
}

impl RegionReport {
    pub fn new() -> Self {
        RegionReport {
            region_name: String::new(),
            filtered_reads_for_region: 0,
            tcs_consensus_results: None,
            umi_summary: None,
        }
    }
}

impl TcsReport {
    pub fn new() -> Self {
        TcsReport {
            process_start_time: Local::now(),
            current_version: env!("CARGO_PKG_VERSION").to_string(),
            input_directory: String::new(),
            advanced_settings: AdvancedSettings::new(),
            input_params: Params::new(),
            total_reads: 0,
            failed_match_reasons: Vec::new(),
            region_reports: Vec::new(),
            errors: Vec::new(),
            warnings: Vec::new(),
            process_end_time: Local::now(),
        }
    }

    pub fn is_successful(&self) -> bool {
        self.errors.is_empty()
    }

    pub fn add_error(&mut self, error: String) {
        self.errors.push(error);
    }

    pub fn add_failed_match_reason(&mut self, reason: FilterPairInvalidReason) {
        self.failed_match_reasons.push(reason);
    }

    pub fn add_warning(&mut self, warning: TcsReportWarnings) {
        self.warnings.push(warning);
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters, Copy)]
pub struct AdvancedSettings {
    #[getset(get = "pub", set = "pub")]
    keep_original: bool,
    #[getset(get = "pub", set = "pub")]
    steepness: f32,
    #[getset(get = "pub", set = "pub")]
    midpoint: u8,
}

impl AdvancedSettings {
    pub fn new() -> Self {
        AdvancedSettings {
            keep_original: false,
            steepness: 0.0,
            midpoint: 0,
        }
    }
    pub fn from_attr(keep_original: bool, steepness: f32, midpoint: u8) -> Self {
        AdvancedSettings {
            keep_original,
            steepness,
            midpoint,
        }
    }

    pub fn default() -> Self {
        AdvancedSettings {
            keep_original: false,
            steepness: 0.2,
            midpoint: 30,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TcsReportWarnings {
    R1R2filteringwarning(String),
    LowAbundanceWarning(String, f64),
    UMIDistErrorWithRegion(String, String),
    ConsensusErrorIndividualWithRegion(String, String),
    EndJoiningErrorWithRegion(String, String),
    QcAndTrimErrorWithRegion(String, String),
}

impl Display for TcsReportWarnings {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            TcsReportWarnings::R1R2filteringwarning(msg) => {
                write!(f, "R1/R2 Filtering Warning: {}", msg)
            }
            TcsReportWarnings::UMIDistErrorWithRegion(region, msg) => {
                write!(f, "UMI distribution error for Region {}: {}", region, msg)
            }
            TcsReportWarnings::ConsensusErrorIndividualWithRegion(region, msg) => {
                write!(
                    f,
                    "Encountered error processing Region: {} for consensus calling, individual consensus aborted, with following error messages: {}",
                    region, msg
                )
            }
            TcsReportWarnings::EndJoiningErrorWithRegion(region, msg) => {
                write!(
                    f,
                    "Encountered error processing Region: {} for end joining, individual end joining aborted, with following error messages: {}",
                    region, msg
                )
            }
            TcsReportWarnings::QcAndTrimErrorWithRegion(region, msg) => {
                write!(
                    f,
                    "Encountered error processing Region: {} for QC and trimming, individual QC and trimming aborted, with following error messages: {}",
                    region, msg
                )
            }
            TcsReportWarnings::LowAbundanceWarning(region, abundance) => {
                write!(
                    f,
                    "Low abundance warning for Region {}: Abundance {} is below threshold: {}, potential cross-contamination",
                    region, abundance, LOW_ABUNDANCE_THRESHOLD_FOR_RAW_READS
                )
            }
        }
    }
}

pub fn tablulate_failed_match_reasons(
    failed_match_reasons: &[FilterPairInvalidReason],
) -> Vec<(FilterPairInvalidReason, usize)> {
    let mut reason_counts: HashMap<FilterPairInvalidReason, usize> = HashMap::new();
    for reason in failed_match_reasons {
        *reason_counts.entry(reason.clone()).or_insert(0) += 1;
    }

    reason_counts
        .into_iter()
        .sorted_by_key(|x| std::cmp::Reverse(x.1))
        .collect::<Vec<_>>()
}

#[derive(Serialize)]
#[serde(untagged)]
pub enum ReasonCount {
    Simple {
        total: usize,
    },
    WithSubcats {
        total: usize,
        #[serde(flatten)]
        subcategories: HashMap<String, usize>,
    },
}

// pub fn tablulate_failed_match_reasons_for_json(
//     failed_match_reasons: &[String],
// ) -> Vec<ReasonCount> {
//     let mut nested_reason_counts = Vec::new();
//     let reason_counts = tablulate_failed_match_reasons(failed_match_reasons);
//     for (main_reason, count) in &reason_counts {
//         if main_reason.contains("No match found") {
//             nested_reason_counts.push(ReasonCount::Simple { total: *count });
//             continue;
//         }

//         let sub = main_reason.split(":").collect::<Vec<&str>>();

//         match sub.as_slice() {
//             [main, subcat] => {
//                 let mut subcategories = HashMap::new();
//                 subcategories.insert(subcat.to_string(), *count);
//                 nested_reason_counts.push(ReasonCount::WithSubcats {
//                     total: *count,
//                     subcategories,
//                 });
//             }
//             [main] => {
//                 nested_reason_counts.push(ReasonCount::Simple { total: *count });
//             }
//             _ => {
//                 // Handle unexpected format
//                 nested_reason_counts.push(ReasonCount::Simple { total: *count });
//             }
//         }
//     }

//     nested_reason_counts
// }
