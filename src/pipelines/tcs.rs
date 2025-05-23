use std::error::Error;
use std::fs::{self, File, OpenOptions};
use std::io::BufWriter;
use std::path::Path;
use std::sync::{Arc, Mutex};

use rayon::prelude::*;

use crate::utils::io::read_fastq_file;
use crate::utils::params::Params;
use crate::utils::tcs_helper::{FastqFiles, log_line, validate_files};

//TODO: write details of the function

pub fn tcs(input: &str, param: &str, keep_original: bool) -> Result<(), Box<dyn Error>> {
    let input_dir = Path::new(input);

    // Check if the input directory exists
    if !input_dir.exists() {
        return Err(format!("Input directory {} does not exist", input).into());
    }
    // Check if the input directory is a valid directory
    if !input_dir.is_dir() {
        return Err(format!("Input path {} is not a valid directory", input).into());
    }

    let logfile = OpenOptions::new()
        .create(true)
        .append(true)
        .open(input_dir.join("run_log.txt"))?;

    let mut logger = BufWriter::new(logfile);

    log_line(&mut logger, "Starting TCS pipeline")?;
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

//TODO: write details of the function
fn run_tcs_pipeline(
    fastq_files: &FastqFiles,
    param: &str,
    logger: &mut BufWriter<File>,
) -> Result<(), Box<dyn Error>> {
    log_line(logger, "Reading Param file")?;
    let params: Params = Params::from_json_sting(&fs::read_to_string(param)?)?;
    log_line(logger, "Validating Params")?;

    // TODO: Validate the params, finish the implementation
    params.validate()?;

    let pairs = read_fastq_file(&fastq_files)?;

    log_line(logger, "Reading Fastq files")?;

    let logs = Arc::new(Mutex::new(Vec::new()));
    let _errors: Arc<Mutex<Vec<Box<dyn std::error::Error + Send + Sync>>>> =
        Arc::new(Mutex::new(Vec::new()));

    // Process the pairs in parallel
    // TODO: currently, this is just a placeholder for the actual processing
    pairs.par_iter().for_each(|(r1, r2)| {
        let id1 = r1.id();
        let id2 = r2.id();
        let msg = if r1.seq().len() == r2.seq().len() {
            format!("Processing pair: {} and {}", id1, id2)
        } else {
            format!("Skipping pair: {} and {}", id1, id2)
        };
        logs.lock().unwrap().push(msg);

        //TODO: populate errors
    });

    for msg in logs.lock().unwrap().iter() {
        log_line(logger, msg)?;
    }

    // TODO: downstream processing
    todo!();
}

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
