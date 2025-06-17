use chrono::{DateTime, Local};
use getset::{Getters, Setters};
use serde::{Deserialize, Serialize};

use crate::helper::tcs_helper::{tcs_consensus::count_joined_and_passed, *};

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters)]
pub struct TcsReportSummary {
    #[getset(get = "pub", set = "pub")]
    process_start_time: DateTime<Local>,
    #[getset(get = "pub", set = "pub")]
    process_end_time: DateTime<Local>,
    #[getset(get = "pub", set = "pub")]
    current_version: String,
    #[getset(get = "pub", set = "pub")]
    input_directory: String,
    #[getset(get = "pub", set = "pub")]
    advanced_settings: AdvancedSettings,
    #[getset(get = "pub", set = "pub")]
    warnings: Vec<TcsReportWarnings>,
    #[getset(get = "pub", set = "pub")]
    region_summaries: Vec<RegionReportSummary>,
}

//TODO: finish this struct
impl TcsReportSummary {
    pub fn new() -> Self {
        TcsReportSummary {
            process_start_time: Local::now(),
            process_end_time: Local::now(),
            current_version: env!("CARGO_PKG_VERSION").to_string(),
            input_directory: String::new(),
            advanced_settings: AdvancedSettings::default(),
            warnings: Vec::new(),
            region_summaries: Vec::new(),
        }
    }

    pub fn from_tcs_report(report: &TcsReport) -> Self {
        let mut summary = TcsReportSummary::new();
        summary.set_process_start_time(*report.process_start_time());
        summary.set_process_end_time(*report.process_end_time());
        summary.set_input_directory(report.input_directory().to_owned());
        summary.set_advanced_settings(report.advanced_settings().clone());
        summary.set_warnings(report.warnings().to_vec());

        for region_report in report.region_reports() {
            let region_summary = RegionReportSummary::from_region_report(region_report);
            summary.region_summaries.push(region_summary);
        }

        summary
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters)]
pub struct RegionReportSummary {
    #[getset(get = "pub", set = "pub")]
    region_name: String,
    #[getset(get = "pub", set = "pub")]
    filtered_reads_for_region: usize,
    #[getset(get = "pub", set = "pub")]
    tcs_number: usize,
    #[getset(get = "pub", set = "pub")]
    umi_cut_off: Option<usize>,
    #[getset(get = "pub", set = "pub")]
    distinct_to_raw_ratio: Option<f64>,
    #[getset(get = "pub", set = "pub")]
    joined_tcs_number: usize,
    #[getset(get = "pub", set = "pub")]
    tcs_passed_qc_number: usize,
}

impl RegionReportSummary {
    pub fn new(region_name: String) -> Self {
        RegionReportSummary {
            region_name,
            filtered_reads_for_region: 0,
            tcs_number: 0,
            distinct_to_raw_ratio: None,
            umi_cut_off: None,
            joined_tcs_number: 0,
            tcs_passed_qc_number: 0,
        }
    }

    //TODO: finish this function
    pub fn from_region_report(region_report: &RegionReport) -> Self {
        let mut region_summary = RegionReportSummary::new(region_report.region_name().to_owned());
        region_summary.set_filtered_reads_for_region(*region_report.filtered_reads_for_region());

        let tcs_consensus_results = region_report.tcs_consensus_results();
        if let Some(results) = tcs_consensus_results {
            region_summary.set_tcs_number(results.len());
            let (n_joined, n_passed) = count_joined_and_passed(results);
            region_summary.set_joined_tcs_number(n_joined);
            region_summary.set_tcs_passed_qc_number(n_passed);
        };

        if let Some(umi_summary) = region_report.umi_summary() {
            region_summary.set_umi_cut_off(Some(*umi_summary.umi_cut_off()));
        }

        if let Some(umi_summary) = region_report.umi_summary() {
            let passed_umis = umi_summary.get_passed_umis_hashmap();
            let distinct_umis = passed_umis.len();
            if distinct_umis > 0 {
                let ratio =
                    distinct_umis as f64 / *region_report.filtered_reads_for_region() as f64;
                region_summary.set_distinct_to_raw_ratio(Some(ratio));
            } else {
                region_summary.set_distinct_to_raw_ratio(Some(0.0));
            }
        } else {
            region_summary.set_distinct_to_raw_ratio(None);
        }

        region_summary
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_tcs_report_summary() {
        let path = "tests/tcs_report.json";
        let report: TcsReport =
            serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
        let summary = TcsReportSummary::from_tcs_report(&report);

        dbg!(&summary);
        assert_eq!(summary.current_version(), env!("CARGO_PKG_VERSION"));
    }
}
