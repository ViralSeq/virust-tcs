use std::collections::HashMap;
use std::error::Error;
use std::fs::File;
use std::io::{BufWriter, Result as IoResult, Write};
use std::path::{self, PathBuf};

use chrono::{DateTime, Local};
use getset::{Getters, Setters};
use serde::{Deserialize, Serialize};

use crate::helper::tcs_helper::*;

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
    total_reads: usize,
    #[getset(get = "pub", set = "pub")]
    warnings: Vec<TcsReportWarnings>,
    #[getset(get = "pub", set = "pub")]
    region_summaries: Vec<RegionReportSummary>,
}

impl TcsReportSummary {
    pub fn new() -> Self {
        TcsReportSummary {
            process_start_time: Local::now(),
            process_end_time: Local::now(),
            current_version: env!("CARGO_PKG_VERSION").to_string(),
            input_directory: String::new(),
            advanced_settings: AdvancedSettings::default(),
            total_reads: 0,
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
        summary.set_total_reads(*report.total_reads());
        summary.set_warnings(report.warnings().to_vec());

        for region_report in report.region_reports() {
            let region_summary = RegionReportSummary::from_region_report(region_report);
            summary.region_summaries.push(region_summary);
        }

        summary
    }

    pub fn to_csv_string(&self) -> Result<String, Box<dyn Error>> {
        let mut wtr = csv::Writer::from_writer(vec![]);

        // Write headers
        wtr.write_record(&[
            "lib_name",
            "region_name",
            "total_reads",
            "filtered_reads_for_region",
            "passed_umis",
            "tcs_number",
            "umi_cut_off",
            "distinct_to_raw_ratio",
            "resampling_index",
            "joined_tcs_number",
            "tcs_passed_qc_number",
        ])?;

        for region in &self.region_summaries {
            wtr.write_record(&[
                PathBuf::from(self.input_directory.clone())
                    .file_name()
                    .unwrap_or_default()
                    .to_string_lossy()
                    .into_owned(),
                region.region_name().to_string(),
                self.total_reads.to_string(),
                region.filtered_reads_for_region().to_string(),
                region.passed_umis().to_string(),
                region.tcs_number().to_string(),
                region
                    .umi_cut_off()
                    .map(|x| x.to_string())
                    .unwrap_or_default(),
                region
                    .distinct_to_raw_ratio()
                    .map(|x| format!("{:.4}", x))
                    .unwrap_or_default(),
                region
                    .resampling_index()
                    .map(|x| format!("{:.4}", x))
                    .unwrap_or_default(),
                region.joined_tcs_number().to_string(),
                region.tcs_passed_qc_number().to_string(),
            ])?;
        }

        wtr.flush()?;
        let data = wtr.into_inner()?; // Vec<u8>
        let csv_string = String::from_utf8(data)?;
        Ok(csv_string)
    }
}

#[derive(Debug, Clone, Serialize, Deserialize, Getters, Setters)]
pub struct RegionReportSummary {
    #[getset(get = "pub", set = "pub")]
    region_name: String,
    #[getset(get = "pub", set = "pub")]
    filtered_reads_for_region: usize,
    #[getset(get = "pub", set = "pub")]
    passed_umis: usize,
    #[getset(get = "pub", set = "pub")]
    tcs_number: usize,
    #[getset(get = "pub", set = "pub")]
    umi_cut_off: Option<usize>,
    #[getset(get = "pub", set = "pub")]
    distinct_to_raw_ratio: Option<f64>,
    #[getset(get = "pub", set = "pub")]
    resampling_index: Option<f64>,
    #[getset(get = "pub", set = "pub")]
    joined_tcs_number: usize,
    #[getset(get = "pub", set = "pub")]
    tcs_passed_qc_number: usize,
    // add a field of detection sensitivity
}

impl RegionReportSummary {
    pub fn new(region_name: String) -> Self {
        RegionReportSummary {
            region_name,
            filtered_reads_for_region: 0,
            passed_umis: 0,
            tcs_number: 0,
            distinct_to_raw_ratio: None,
            umi_cut_off: None,
            resampling_index: None,
            joined_tcs_number: 0,
            tcs_passed_qc_number: 0,
        }
    }

    pub fn from_region_report(region_report: &RegionReport) -> Self {
        let mut region_summary = RegionReportSummary::new(region_report.region_name().to_owned());
        region_summary.set_filtered_reads_for_region(*region_report.filtered_reads_for_region());

        let tcs_consensus_results = region_report.tcs_consensus_results();
        if let Some(results) = tcs_consensus_results {
            region_summary.set_tcs_number(results.len());
            let (n_joined, n_passed) = tcs_consensus::count_joined_and_passed(results);
            region_summary.set_joined_tcs_number(n_joined);
            region_summary.set_tcs_passed_qc_number(n_passed);
        };

        if let Some(umi_summary) = region_report.umi_summary() {
            region_summary.set_umi_cut_off(Some(*umi_summary.umi_cut_off()));
        }

        if let Some(umi_summary) = region_report.umi_summary() {
            let passed_umis = umi_summary.get_passed_umis_hashmap();
            region_summary.set_passed_umis(passed_umis.len());

            let distinct_umis = passed_umis.len();
            if distinct_umis > 0 {
                let ratio =
                    distinct_umis as f64 / *region_report.filtered_reads_for_region() as f64;
                region_summary.set_distinct_to_raw_ratio(Some(ratio));

                if passed_umis.len() > 0 {
                    region_summary.set_resampling_index(Some(
                        *region_summary.tcs_number() as f64 / passed_umis.len() as f64,
                    ));
                } else {
                    region_summary.set_resampling_index(None);
                }
            } else {
                region_summary.set_distinct_to_raw_ratio(Some(0.0));
            }
        } else {
            region_summary.set_distinct_to_raw_ratio(None);
        }

        region_summary
    }
}

pub fn tcs_report_write(
    tcs_report: &TcsReport,
    report_logger: &mut BufWriter<File>,
) -> IoResult<()> {
    writeln!(
        report_logger,
        "{}",
        serde_json::to_string_pretty(&TcsReportSummary::from_tcs_report(&tcs_report))? // if the full raw report is needed, use the following line instead:
                                                                                       // serde_json::to_string_pretty(&tcs_report)?
    )?;
    report_logger.flush()?;

    Ok(())
}

pub fn raw_sequence_invalid_reason_write(
    tcs_report: &TcsReport,
    path: &str,
) -> Result<(), Box<dyn Error>> {
    let output_path = path::Path::new(path);
    if !output_path.exists() {
        return Err(
            TcsError::UnexpectedError("Unable to access the output directory".to_string()).into(),
        );
    }

    let file = output_path.join("raw_sequence_invalid_reasons.csv");

    let reasons = tablulate_failed_match_reasons(&tcs_report.failed_match_reasons());
    let data = serde_json::to_string_pretty(&reasons)?;
    let json: Vec<(HashMap<String, String>, usize)> = serde_json::from_str(&data)?;

    // 1. first, sum totals by main category
    let mut main_totals: HashMap<String, usize> = HashMap::new();
    for (cat_map, count) in &json {
        let main_cat = cat_map.keys().next().unwrap();
        *main_totals.entry(main_cat.clone()).or_insert(0) += count;
    }

    // 2. Now, output rows with all four columns

    let mut writer = csv::Writer::from_path(file)?;
    writer.write_record(&["main_category", "sub_category", "count_sub", "count_main"])?;

    for (cat_map, count) in &json {
        let main_cat = cat_map.keys().next().unwrap();
        let sub_cat = cat_map.values().next().unwrap();
        let main_count = main_totals.get(main_cat).unwrap_or(&0);
        writer.write_record(&[
            main_cat,
            sub_cat,
            &count.to_string(),
            &main_count.to_string(),
        ])?;
    }
    writer.flush()?;
    Ok(())
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
        assert_eq!(summary.region_summaries().len(), 5);
    }
}
