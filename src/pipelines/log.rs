//TODO Log pipeline
use std::error::Error;
use std::fs;
use std::path::PathBuf;

use crate::helper::io::find_directories;
use crate::helper::json::FromJsonString;
use crate::helper::tcs_helper::*;

pub fn run_log(input: String, output: String) -> Result<(), Box<dyn Error>> {
    let directories = find_directories(&input)?;

    let lib_names = directories
        .iter()
        .map(|path| path.file_name().unwrap().to_string_lossy())
        .collect::<Vec<_>>();
    println!("Found directories: {:?}", lib_names);

    let mut summaries: Vec<TcsReportSummary> = Vec::new();

    for dir in directories {
        let summary_file_path = dir.join("tcs_report.json");
        if !summary_file_path.exists() {
            println!("No TCS summary file found in directory: {}", dir.display());
            continue;
        }

        let tcs_summary = tcs_summary::TcsReportSummary::from_json_string(
            &std::fs::read_to_string(&summary_file_path)?,
        )?;

        summaries.push(tcs_summary);
    }
    let final_tcs_summary_csv_report = merge_csv_summaries(&summaries)?;

    let output_path = PathBuf::from(output);

    if output_path.is_file() {
        return Err("Output path must be a directory".to_string().into());
    } else if !output_path.exists() {
        std::fs::create_dir_all(&output_path)?;
    }

    let csv_log_file = output_path.join("log.csv");

    fs::write(csv_log_file, final_tcs_summary_csv_report)?;

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
