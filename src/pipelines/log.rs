//TODO Log pipeline
use std::error::Error;
use std::fs;
use std::path::PathBuf;

use bio::io::fastq::Record;
use bio::io::{fasta, fastq};
use flate2::Compression;
use flate2::write::GzEncoder;

use crate::helper::fastqc;
use crate::helper::io::find_directories;
use crate::helper::json::FromJsonString;
use crate::helper::params::Params;
use crate::helper::tcs_helper::*;
use crate::helper::umis;

pub fn run_log(input: String, output: String) -> Result<(), Box<dyn Error>> {
    let output_path = PathBuf::from(output);

    if output_path.is_file() {
        return Err("Output path must be a directory".to_string().into());
    } else if !output_path.exists() {
        fs::create_dir_all(&output_path)?;
    }

    let fasta_dir = output_path.join("Joined_TCS_fasta");
    let fastq_dir = output_path.join("Joined_TCS_fastq");
    let fastq_qc_dir = output_path.join("FastQC_reports");
    let temp_data_dir = output_path.join("temp_data");

    if !fasta_dir.exists() {
        fs::create_dir_all(&fasta_dir)?;
    }
    if !fastq_dir.exists() {
        fs::create_dir_all(&fastq_dir)?;
    }
    if !fastq_qc_dir.exists() {
        fs::create_dir_all(&fastq_qc_dir)?;
    }
    if !temp_data_dir.exists() {
        fs::create_dir_all(&temp_data_dir)?;
    }

    let directories = find_directories(&input)?;

    let mut summaries: Vec<TcsReportSummary> = Vec::new();

    for dir in directories {
        let lib_name = dir.file_name().unwrap().to_string_lossy();
        println!("Processing directory: {} ({})", lib_name, dir.display());
        let fasta_dir_with_lib = fasta_dir.join(lib_name.as_ref());
        let fastq_dir_with_lib = fastq_dir.join(lib_name.as_ref());
        let fastq_qc_dir_with_lib = fastq_qc_dir.join(lib_name.as_ref());
        let temp_data_dir_with_lib = temp_data_dir.join(lib_name.as_ref());
        if !fasta_dir_with_lib.exists() {
            fs::create_dir_all(&fasta_dir_with_lib)?;
        }
        if !fastq_dir_with_lib.exists() {
            fs::create_dir_all(&fastq_dir_with_lib)?;
        }
        let summary_file_path = dir.join("tcs_report.json");
        if !summary_file_path.exists() {
            println!("No TCS summary file found in directory: {}", dir.display());
            continue;
        }
        if !fastq_qc_dir_with_lib.exists() {
            fs::create_dir_all(&fastq_qc_dir_with_lib)?;
        }
        if !temp_data_dir_with_lib.exists() {
            fs::create_dir_all(&temp_data_dir_with_lib)?;
        }

        let tcs_summary = tcs_summary::TcsReportSummary::from_json_string(&fs::read_to_string(
            &summary_file_path,
        )?)?;

        summaries.push(tcs_summary.clone());

        let sample_log_file = temp_data_dir_with_lib.join("sample_log.csv");

        fs::write(sample_log_file, merge_csv_summaries(&[tcs_summary])?)?;

        let params = Params::from_json_string(&fs::read_to_string(dir.join("tcs_params.json"))?)?;

        let rsfr_file = dir.join("raw_sequence_invalid_reasons.csv");
        if rsfr_file.exists() {
            fs::copy(
                &rsfr_file,
                temp_data_dir_with_lib.join("raw_sequence_invalid_reasons.csv"),
            )?;
        }

        let umi_dis_file = temp_data_dir_with_lib.join("umi_distribution.csv");

        let mut umi_wtr = csv::Writer::from_path(umi_dis_file)?;

        umi_wtr.write_record(&["region", "umi", "umi_count", "umi_cut_off"])?;

        // get directorys from this path
        let subdirectories = find_directories(dir.to_str().unwrap())?;
        for subdir in subdirectories {
            let region_name = subdir.file_name().unwrap().to_string_lossy();
            let joined_fastq_name = determine_joined_tcs_file_from_params(&params, &region_name);
            if joined_fastq_name.is_none() {
                println!("Region: {}, No joined FASTQ found", region_name);
                continue;
            }
            let joined_fastq_name = find_fastq(&subdir, &joined_fastq_name.unwrap());

            if let Some(joined_fastq) = joined_fastq_name {
                let fastq_reader = fastq::Reader::from_file(&joined_fastq)?;
                let mut fasta_writer = fasta::Writer::to_file(
                    fasta_dir_with_lib.join(format!("{}_{}.fasta", lib_name, region_name)),
                )?;
                let mut fastq_writer = fastq::Writer::to_file(
                    fastq_dir_with_lib.join(format!("{}_{}.fastq", lib_name, region_name)),
                )?;

                for record in fastq_reader.records() {
                    let record = record?;
                    let new_id = format!("{}|{}|{}", lib_name, region_name, record.id());
                    let new_record =
                        Record::with_attrs(&new_id, record.desc(), record.seq(), record.qual());
                    fastq_writer.write_record(&new_record)?;
                    fasta_writer.write_record(&fasta::Record::with_attrs(
                        &new_id,
                        record.desc(),
                        record.seq(),
                    ))?;
                }

                drop(fastq_writer);
                drop(fasta_writer);

                // run fastqc analysis
                let fastqc_results = fastqc::fastqc_analysis(&joined_fastq)?;
                let qc_report_path =
                    fastq_qc_dir_with_lib.join(format!("{}_{}_fastqc.png", lib_name, region_name));
                fastqc::plot_quality_score_distribution(
                    &fastqc_results.quality_score_distribution(),
                    &qc_report_path,
                )?;
                fastqc_results.export_quality_score_distribution_to_csv(
                    &fastq_qc_dir_with_lib.join(format!("{}_{}_fastqc.csv", lib_name, region_name)),
                )?;

                // compress the joined fastq, and remove the original uncompressed file
                compress_fastq_gz(
                    &fastq_dir_with_lib.join(format!("{}_{}.fastq", lib_name, region_name)),
                )?;

                let umi_summary_file = subdir.join("umi_summary.json");
                if umi_summary_file.exists() {
                    let umi_summary = umis::UMISummary::from_json_string(&fs::read_to_string(
                        &umi_summary_file,
                    )?)?;

                    for (umi, umi_count) in umi_summary.umi_freq() {
                        umi_wtr.write_record(vec![
                            region_name.to_string(),
                            umi.to_string(),
                            umi_count.to_string(),
                            umi_summary.umi_cut_off().to_string(),
                        ])?;
                    }
                }
            } else {
                println!("Region: {}, No joined FASTQ found", region_name);
            }
        }

        umi_wtr.flush()?;
    }

    let final_tcs_summary_csv_report = merge_csv_summaries(&summaries)?;

    let csv_log_file = output_path.join("log.csv");

    fs::write(&csv_log_file, final_tcs_summary_csv_report)?;

    fs::copy(&csv_log_file, temp_data_dir.join("log.csv"))?;

    Ok(())
}

// Merges multiple TCS report summaries into a single CSV string
fn merge_csv_summaries(summaries: &[TcsReportSummary]) -> Result<String, Box<dyn Error>> {
    let mut final_csv_report = String::new();

    for (i, summary) in summaries.iter().enumerate() {
        let csv_str = summary.to_csv_string()?;

        if i == 0 {
            // first one: keep hearder
            final_csv_report.push_str(&csv_str);
        } else {
            // skip the header line for subsequent
            if let Some(pos) = csv_str.find('\n') {
                final_csv_report.push_str(&csv_str[pos + 1..]);
            }
        }
    }
    Ok(final_csv_report)
}

// Find the matching FASTQ under the fastq_files/ within a TCS/Region output directory
fn find_fastq(root: &PathBuf, target_name: &str) -> Option<PathBuf> {
    let candidate = root.join("fastq_files").join(target_name);
    if candidate.exists() {
        Some(candidate)
    } else {
        None
    }
}

fn determine_joined_tcs_file_from_params(params: &Params, region_name: &str) -> Option<String> {
    for region_param in &params.primer_pairs {
        if region_param.region == region_name {
            if region_param.trim {
                return Some("joined_passed_qc_trimmed.fastq".to_string());
            } else if region_param.tcs_qc {
                return Some("joined_passed_qc.fastq".to_string());
            } else {
                return Some("joined.fastq".to_string());
            }
        }
    }
    None
}

fn compress_fastq_gz(input: &PathBuf) -> Result<(), Box<dyn Error>> {
    let mut output = input.clone();
    output.set_extension("fastq.gz");
    let input_file = fs::File::open(&input)?;
    let output_file = fs::File::create(&output)?;
    let mut encoder = GzEncoder::new(output_file, Compression::default());
    std::io::copy(&mut std::io::BufReader::new(input_file), &mut encoder)?;
    encoder.finish()?;

    fs::remove_file(input)?; // Remove the original uncompressed file
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_determine_joined_tcs_file_from_params() {
        let path_to_params = "tests/data/pr.json";
        let params =
            Params::from_json_string(&std::fs::read_to_string(path_to_params).unwrap()).unwrap();
        let region_name = "PR";
        let result = determine_joined_tcs_file_from_params(&params, region_name);

        assert!(result.is_some());
        assert_eq!(result.unwrap(), "joined_passed_qc_trimmed.fastq");
    }
}
