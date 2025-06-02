use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File, OpenOptions};
use std::io::BufWriter;
use std::path::Path;

use rayon::prelude::*;

use crate::helper::consensus::*;
use crate::helper::io::read_fastq_file;
use crate::helper::params::Params;
use crate::helper::tcs_helper::TcsError;
use crate::helper::tcs_helper::log_line;
use crate::helper::tcs_helper::{
    FastqFiles, FilteredPair, PairedRecordFilterResult, filter_r1_r2_pairs, validate_files,
};

// MARK: tcs main function
//TODO: write details of the function

pub fn tcs(
    input: &str,
    param: &str,
    keep_original: bool,
    steepness: f32,
    midpoint: u8,
) -> Result<(), Box<dyn Error>> {
    let input_dir = Path::new(input);

    // Check if the input directory exists
    if !input_dir.exists() {
        return Err(TcsError::InputDirNotFound(input.to_string()).into());
    }
    // Check if the input directory is a valid directory
    if !input_dir.is_dir() {
        return Err(TcsError::NotADirectory(input.to_string()).into());
    }

    let logfile = OpenOptions::new()
        .create(true)
        .append(true)
        .open(input_dir.join("run_log.txt"))?;

    let mut logger = BufWriter::new(logfile);

    log_line(&mut logger, "Starting TCS pipeline")?;
    log_line(
        &mut logger,
        &format!("TCS (Rust) Version: {}", env!("CARGO_PKG_VERSION")),
    )?;
    log_line(&mut logger, &format!("Input directory: {}", input))?;
    log_line(&mut logger, &format!("Param file input: {}", param))?;
    log_line(&mut logger, &format!("Keep original: {}", keep_original))?;
    log_line(&mut logger, &format!("Steepness: {}", steepness))?;
    log_line(&mut logger, &format!("Midpoint: {}", midpoint))?;
    log_line(&mut logger, "Validating input files")?;

    let fastq_files = validate_files(input)?;
    let r1_file = &fastq_files.r1_file;
    let r2_file = &fastq_files.r2_file;
    let data_type = &fastq_files.data_type;
    log_line(&mut logger, &format!("R1 file: {:?}", r1_file))?;
    log_line(&mut logger, &format!("R2 file: {:?}", r2_file))?;
    log_line(&mut logger, &format!("Data type: {:?}", data_type))?;

    if let Err(e) = run_tcs_pipeline(&fastq_files, param, &mut logger, steepness, midpoint) {
        log_line(&mut logger, &format!("Error running TCS pipeline: {}", e))?;
        return Err(e);
    }

    if keep_original {
        log_line(&mut logger, "Keeping original files")?;
    } else {
        log_line(&mut logger, "Deleting original files")?;
        std::fs::remove_file(r1_file)?;
        std::fs::remove_file(r2_file)?;
    }
    log_line(&mut logger, "TCS pipeline completed")?;

    Ok(())
}

// MARK: run_tcs_pipeline function
//TODO: write details of the function
fn run_tcs_pipeline(
    fastq_files: &FastqFiles,
    param: &str,
    logger: &mut BufWriter<File>,
    steepness: f32,
    midpoint: u8,
) -> Result<(), Box<dyn Error>> {
    log_line(logger, "Reading Param file")?;
    let params: Params = Params::from_json_sting(&fs::read_to_string(param)?)?;
    log_line(logger, "Validating Params")?;

    let validated_params = params.validate()?;

    let pairs = read_fastq_file(&fastq_files)?;

    log_line(logger, "Reading Fastq files")?;
    log_line(
        logger,
        &format!("Number of raw fastq records: {}", pairs.len()),
    )?;

    // Process the pairs in parallel

    let (groups, fails, errors) = pairs
        .par_iter()
        .fold(
            // Each thread starts with its own empty results
            || (HashMap::new(), Vec::new(), Vec::new()),
            |(mut ok, mut fail, mut err), pair| {
                match filter_r1_r2_pairs(&pair.0, &pair.1, &validated_params) {
                    Ok(filter_result) => match filter_result {
                        PairedRecordFilterResult::Valid(filtered_pair) => {
                            let region = filtered_pair.region.clone();
                            ok.entry(region)
                                .or_insert_with(Vec::new)
                                .push(filtered_pair);
                        }
                        PairedRecordFilterResult::Invalid(reason) => {
                            fail.push(reason);
                        }
                    },
                    Err(e) => {
                        err.push(e);
                    }
                }
                (ok, fail, err)
            },
        )
        .reduce(
            // Combine the results from all threads
            || {
                (
                    HashMap::<String, Vec<FilteredPair>>::new(),
                    Vec::new(),
                    Vec::new(),
                )
            },
            |(mut ok1, mut fail1, mut err1), (ok2, fail2, err2)| {
                for (region, mut vec) in ok2 {
                    ok1.entry(region).or_insert_with(Vec::new).append(&mut vec);
                }
                fail1.extend(fail2);
                err1.extend(err2);
                (ok1, fail1, err1)
            },
        );

    let mut fail_frequency = HashMap::new();
    for fail in &fails {
        *fail_frequency.entry(fail.clone()).or_insert(0) += 1;
    }

    for fail_freq in fail_frequency {
        log_line(
            logger,
            &format!("Fail reason: {}, count: {}", fail_freq.0, fail_freq.1),
        )?;
    }
    for error in errors {
        log_line(logger, &format!("Error: {}", error))?;
    }

    for (region, filtered_pairs) in groups {
        log_line(
            logger,
            &format!(
                "Region: {}, valid r1 r2 pairs: {}",
                region,
                filtered_pairs.len()
            ),
        )?;
    }

    let consensus_params = ConsensusParams::new(steepness as f64, midpoint as f64);

    // TODO: downstream processing
    // 1. consensus calling for each region.
    // 2. end-joining of consensus sequences.
    // 3. qc
    // 4. trimming
    // 5. write fastq and fasta files.
    // 6. write UMI files in JSON format.
    // 7. export a summary report.
    todo!(
        "Downstream processing not implemented yet. Consensus params: k {}, q0 {}",
        consensus_params.k(),
        consensus_params.q0()
    );
}

// MARK: tcs_dr main function
//TODO: write details of the function
pub fn tcs_dr(
    input: &str,
    version: Option<String>,
    keep_original: bool,
) -> Result<(), Box<dyn Error>> {
    todo!(
        "TCS DR function called with input: {}, version: {:?}, keep_original: {}",
        input,
        version,
        keep_original
    );
}
