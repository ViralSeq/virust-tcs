use super::*;

use std::error::Error;
use std::fs;
use std::fs::File;
use std::ops::Range;
use std::path::Path;

use bio::io::{fasta, fastq};
use getset::{Getters, Setters};

#[derive(Debug, Clone, Getters, Setters)]
pub struct TcsOutput<'a> {
    #[getset(get = "pub")]
    r1_fastq: Vec<&'a fastq::Record>,
    #[getset(get = "pub")]
    r2_fastq: Vec<&'a fastq::Record>,
    #[getset(get = "pub")]
    joined_tcs_fastq: Option<Vec<&'a fastq::Record>>,
    #[getset(get = "pub")]
    joined_tcs_passed_qc_fastq: Option<Vec<&'a fastq::Record>>,
    #[getset(get = "pub")]
    qc_failed_reasons: Option<Vec<(String, QcNotPassedReport)>>,
}

impl<'a> TcsOutput<'a> {
    pub fn from_region_report(region_report: &'a RegionReport) -> Self {
        if region_report.tcs_consensus_results().is_none() {
            return TcsOutput {
                r1_fastq: Vec::new(),
                r2_fastq: Vec::new(),
                joined_tcs_fastq: None,
                joined_tcs_passed_qc_fastq: None,
                qc_failed_reasons: None,
            };
        }

        let mut r1_seqs = Vec::new();
        let mut r2_seqs = Vec::new();
        let mut joined_seqs = Vec::new();
        let mut joined_passed_qc_seqs = Vec::new();
        let mut qc_failed_reasons = Vec::new();

        for tcs in region_report.tcs_consensus_results().as_ref().unwrap() {
            r1_seqs.push(tcs.r1_consensus());
            r2_seqs.push(tcs.r2_consensus());

            if let Some(joined) = tcs.joined_consensus() {
                joined_seqs.push(joined);
            }

            if *tcs.qc() == TcsConsensusQcResult::Passed {
                if let Some(joined) = tcs.joined_consensus() {
                    joined_passed_qc_seqs.push(joined);
                }
            }
            if let TcsConsensusQcResult::NotPassed(reason) = tcs.qc() {
                qc_failed_reasons.push((tcs.umi_information_block().to_string(), reason.clone()));
            }
        }

        TcsOutput {
            r1_fastq: r1_seqs,
            r2_fastq: r2_seqs,
            joined_tcs_fastq: if joined_seqs.is_empty() {
                None
            } else {
                Some(joined_seqs)
            },
            joined_tcs_passed_qc_fastq: if joined_passed_qc_seqs.is_empty() {
                None
            } else {
                Some(joined_passed_qc_seqs)
            },
            qc_failed_reasons: if qc_failed_reasons.is_empty() {
                None
            } else {
                Some(qc_failed_reasons)
            },
        }
    }
}

pub fn tcs_sequence_data_write(tcs_report: &TcsReport, path: &str) -> Result<(), Box<dyn Error>> {
    let output_path = Path::new(path);
    if !output_path.exists() {
        return Err(
            TcsError::UnexpectedError("Unable to access the output directory".to_string()).into(),
        );
    }

    for region_report in tcs_report.region_reports() {
        let region_dir = output_path.join(region_report.region_name());
        fs::create_dir_all(region_dir.join("fastq_files"))?;
        fs::create_dir_all(region_dir.join("fasta_files"))?;

        let umi_summary_file = region_dir.join("umi_summary.json");

        // write UMI summary if it exists
        if let Some(umi_summary) = region_report.umi_summary() {
            let umi_summary_json = serde_json::to_string_pretty(umi_summary)?;
            fs::write(umi_summary_file, umi_summary_json)?;
        } else {
            // If no UMI summary, create an empty file
            fs::write(umi_summary_file, "{}")?;
        }

        let tcs_output = TcsOutput::from_region_report(region_report);

        write_fastq_and_fasta(&tcs_output.r1_fastq(), &region_dir, "r1")?;
        write_fastq_and_fasta(&tcs_output.r2_fastq(), &region_dir, "r2")?;

        if let Some(joined) = tcs_output.joined_tcs_fastq() {
            write_fastq_and_fasta(joined, &region_dir, "joined")?;
        }
        if let Some(joined_passed_qc) = tcs_output.joined_tcs_passed_qc_fastq() {
            write_fastq_and_fasta(joined_passed_qc, &region_dir, "joined_passed_qc")?;
        }

        let qc_failed_reasons_file = region_dir.join("qc_failed_reasons.csv");
        let mut csv_writer = csv::Writer::from_path(qc_failed_reasons_file)?;
        csv_writer.write_record([
            "UMI",
            "qc_reference",
            "qc_coordinates1_start",
            "qc_coordinates1_end",
            "qc_coordinates2_start",
            "qc_coordinates2_end",
            "qc_indels",
            "locator_coordinates_start",
            "locator_coordinates_end",
            "locator_indels",
        ])?;

        if let Some(qc_reasons) = tcs_output.qc_failed_reasons() {
            for (umi, reason) in qc_reasons {
                let (coord1_start, coord1_end) = match_coordinates(reason.qc_coordinates1());
                let (coord2_start, coord2_end) = match_coordinates(reason.qc_coordinates2());
                let (loc_start, loc_end) = match_coordinates(reason.locator_coordinates());
                csv_writer.write_record([
                    umi,
                    &reason.qc_reference().to_string(),
                    &coord1_start.to_string(),
                    &coord1_end.to_string(),
                    &coord2_start.to_string(),
                    &coord2_end.to_string(),
                    &reason.qc_indels().to_string(),
                    &loc_start.to_string(),
                    &loc_end.to_string(),
                    &reason.locator_indels().to_string(),
                ])?;
            }
        }

        csv_writer.flush()?;
    }
    Ok(())
}

pub fn export_input_params(
    tcs_report: &TcsReport,
    input_directory: &str,
) -> Result<(), Box<dyn Error>> {
    let params = tcs_report.input_params();
    let params_file_path = Path::new(input_directory).join("tcs_params.json");
    let params_file = File::create(params_file_path)?;
    serde_json::to_writer_pretty(params_file, &params)?;
    Ok(())
}

fn write_fastq_and_fasta(
    records: &[&fastq::Record],
    region_dir: &Path,
    prefix: &str,
) -> Result<(), Box<dyn Error>> {
    let fastq_path = region_dir
        .join("fastq_files")
        .join(format!("{prefix}.fastq"));
    let fasta_path = region_dir
        .join("fasta_files")
        .join(format!("{prefix}.fasta"));

    let mut fastq_writer = fastq::Writer::to_file(&fastq_path)?;
    let mut fasta_writer = fasta::Writer::to_file(&fasta_path)?;

    for record in records {
        fastq_writer.write_record(record)?;
        fasta_writer.write_record(&fastq_to_fasta_record(record))?;
    }
    Ok(())
}

fn match_coordinates(coord: &Option<Range<u32>>) -> (u32, u32) {
    match coord {
        Some(range) => (range.start, range.end),
        None => (0, 0), // Default values if no coordinates are provided
    }
}
