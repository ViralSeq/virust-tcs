use std::fmt::Display;

use chrono::{DateTime, Local};
use getset::{Getters, Setters};
use serde::{Deserialize, Serialize};

use crate::helper::params::Params;
use crate::helper::tcs_helper::TcsConsensus;
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
    failed_match_reasons: Vec<String>,

    #[getset(get = "pub", set = "pub")]
    region_reports: Vec<RegionReport>,

    errors: Vec<String>,
    warnings: Vec<TcsReportWarnings>,
    process_end_time: DateTime<Local>,
}

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters)]
pub struct RegionReport {
    #[getset(get = "pub", set = "pub")]
    region_name: String,
    #[getset(get = "pub", set = "pub")]
    filtered_reads_for_region: usize,
    #[getset(get = "pub", set = "pub")]
    tcs_consensus_results: Option<TcsConsensus>,
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

    pub fn add_failed_match_reason(&mut self, reason: String) {
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
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum TcsReportWarnings {
    R1R2filteringwarning(String),
    UMIDistErrorWithRegion(String, String),
    ConsensusErrorIndividualWithRegion(String, String),
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
        }
    }
}
