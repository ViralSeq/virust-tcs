use std::error::Error;
use std::fs::{self, File, OpenOptions};
use std::io::BufWriter;
use std::path::Path;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;

use crate::helper::io::read_fastq_file;
use crate::helper::params::Params;
use crate::helper::tcs_helper::TcsError;
use crate::helper::tcs_helper::log_line;
use crate::helper::tcs_helper::{
    FastqFiles, PairedRecordFilterResult, filter_r1_r2_pairs, validate_files,
};

// MARK: tcs main function
//TODO: write details of the function

pub fn tcs(input: &str, param: &str, keep_original: bool) -> Result<(), Box<dyn Error>> {
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
    log_line(&mut logger, "Validating input files")?;

    let fastq_files = validate_files(input)?;
    let r1_file = &fastq_files.r1_file;
    let r2_file = &fastq_files.r2_file;
    let data_type = &fastq_files.data_type;
    log_line(&mut logger, &format!("R1 file: {:?}", r1_file))?;
    log_line(&mut logger, &format!("R2 file: {:?}", r2_file))?;
    log_line(&mut logger, &format!("Data type: {:?}", data_type))?;

    if let Err(e) = run_tcs_pipeline(&fastq_files, param, &mut logger) {
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
) -> Result<(), Box<dyn Error>> {
    log_line(logger, "Reading Param file")?;
    let params: Params = Params::from_json_sting(&fs::read_to_string(param)?)?;
    log_line(logger, "Validating Params")?;

    let validated_params = params.validate()?;

    let pairs = read_fastq_file(&fastq_files)?;

    log_line(logger, "Reading Fastq files")?;

    let logs = Arc::new(Mutex::new(Vec::new()));

    let errors: Arc<Mutex<Vec<Box<dyn Error + Send + Sync>>>> = Arc::new(Mutex::new(Vec::new()));
    // Process the pairs in parallel
    // TODO: currently, this is just a placeholder for the actual processing
    pairs.par_iter().for_each(|(r1, r2)| {
        match filter_r1_r2_pairs(r1, r2, &validated_params) {
            Ok(filtered_pair) => {
                if let PairedRecordFilterResult::Valid(filtered_pair) = filtered_pair {
                    // Log the filtered pair
                    let msg = format!(
                        "Filtered pair: R1: {}, R2: {}",
                        filtered_pair.r1.id(),
                        filtered_pair.r1.id()
                    );
                    logs.lock().unwrap().push(msg);
                } else {
                    logs.lock()
                        .unwrap()
                        .push("Filtered pair is None".to_string());
                }
            }
            Err(e) => {
                // Log the error
                let mut errs = errors.lock().unwrap();
                errs.push(e);
            }
        }
    });

    for msg in logs.lock().unwrap().iter() {
        log_line(logger, msg)?;
    }

    // TODO: downstream processing
    todo!();
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
