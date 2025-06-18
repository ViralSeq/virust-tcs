use super::*;

use std::error::Error;
use std::fs;
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
}

impl<'a> TcsOutput<'a> {
    pub fn from_region_report(region_report: &'a RegionReport) -> Self {
        if region_report.tcs_consensus_results().is_none() {
            return TcsOutput {
                r1_fastq: Vec::new(),
                r2_fastq: Vec::new(),
                joined_tcs_fastq: None,
                joined_tcs_passed_qc_fastq: None,
            };
        }

        let mut r1_seqs = Vec::new();
        let mut r2_seqs = Vec::new();
        let mut joined_seqs = Vec::new();
        let mut joined_passed_qc_seqs = Vec::new();

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
    }
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
