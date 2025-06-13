use std::collections::HashMap;
use std::error::Error;
use std::fs::{self, File, OpenOptions};
use std::io::BufWriter;
use std::path::{Path, PathBuf};

use rayon::prelude::*;

use crate::helper::consensus::*;
use crate::helper::io::read_fastq_file;
use crate::helper::params::Params;
use crate::helper::tcs_helper::*;

pub fn tcs(
    input: &str,
    param: &str,
    keep_original: bool,
    steepness: f32,
    midpoint: u8,
) -> Result<(), Box<dyn Error>> {
    // initialize the TCS report and logger
    // this will create a new TCS report and a logger that will log the progress of the TCS pipeline.
    // the logger will log to a file named run_log.txt in the input directory.
    // the TCS report will be used to store the results of the TCS pipeline.
    let (mut tcs_report, mut logger) = tcs_init(input)?;

    let advanced_settings = AdvancedSettings::from_attr(keep_original, steepness, midpoint);
    tcs_report.set_advanced_settings(advanced_settings);

    // log the start of the TCS pipeline
    log_line(&mut logger, "Starting TCS pipeline")?;

    // Run the TCS main function

    let (tcs_report, r1_r2_path) = match tcs_main(tcs_report, &mut logger, param, advanced_settings)
    {
        Ok((report, r1_r2_path)) => {
            log_line(&mut logger, "TCS main function completed successfully")?;
            (report, r1_r2_path)
        }
        Err(e) => {
            log_line(
                &mut logger,
                &format!("Fatal error in TCS main function: {}", e),
            )?;
            return Err(e);
        }
    };

    let success = tcs_report.is_successful();

    dbg!(tcs_report.is_successful());

    // TODO: write the TCS report to a file
    // tcs_write(&tcs_report, &mut logger)?;

    if keep_original {
        log_line(&mut logger, "Keeping original files")?;
    } else if r1_r2_path.is_some() {
        log_line(&mut logger, "Deleting original files")?;
        let r1_r2_path = r1_r2_path.unwrap();
        std::fs::remove_file(r1_r2_path.0)?;
        std::fs::remove_file(r1_r2_path.1)?;
    } else {
        log_line(&mut logger, "No original files to delete")?;
    }

    if success {
        log_line(&mut logger, "TCS pipeline completed successfully\n")?;
    } else {
        log_line(&mut logger, "TCS pipeline completed with errors\n")?;
    }

    Ok(())
}

// MARK: tcs main function
//TODO: write details of the function

pub fn tcs_main(
    mut tcs_report: TcsReport,
    logger: &mut BufWriter<File>,
    param: &str,
    advanced_settings: AdvancedSettings,
) -> Result<(TcsReport, Option<(PathBuf, PathBuf)>), Box<dyn Error>> {
    let keep_original = *advanced_settings.keep_original();
    let steepness = *advanced_settings.steepness();
    let midpoint = *advanced_settings.midpoint();
    let input = tcs_report.input_directory();

    log_line(
        logger,
        &format!("TCS (Rust) Version: {}", env!("CARGO_PKG_VERSION")),
    )?;
    log_line(logger, &format!("Input directory: {}", input))?;
    log_line(logger, &format!("Param file input: {}", param))?;
    log_line(logger, &format!("Keep original: {}", keep_original))?;
    log_line(logger, &format!("Steepness: {}", steepness))?;
    log_line(logger, &format!("Midpoint: {}", midpoint))?;
    log_line(logger, "Validating input files")?;

    // Validate the input files and get the fastq files
    // This will check if the input files are valid and return a FastqFiles struct containing the paths to the R1 and R2 files.
    // If the input files are not valid, it will return an error.
    // The FastqFiles struct will also contain the data type (Fastq or FastqGz) of the files.
    // The function validate_files will also log the input files and data type to the logger.
    // If there is an error validating the input files, it will log the error to the logger and return a TcsReport with the error.
    // The TcsReport with error will be handled in the downstream processing.
    let fastq_files = match validate_files(input) {
        Ok(files) => files,
        Err(e) => {
            log_line(logger, &format!("Error validating input files: {}", e))?;
            tcs_report.add_error(e.to_string());
            return Ok((tcs_report, None));
        }
    };

    let r1_file = &fastq_files.r1_file;
    let r2_file = &fastq_files.r2_file;
    let data_type = &fastq_files.data_type;
    log_line(logger, &format!("R1 file: {:?}", r1_file))?;
    log_line(logger, &format!("R2 file: {:?}", r2_file))?;
    log_line(logger, &format!("Data type: {:?}", data_type))?;

    // Read the param file and validate it
    // This will read the param file and parse it into a Params struct.
    // If the param file is not valid, it will return an error.
    // The Params struct will contain the parameters for the TCS pipeline.
    // If there is an error reading the param file, it will log the error to the logger and return a TcsReport with the error.
    // The TcsReport with error will be handled in the downstream processing.
    log_line(logger, "Reading Param file")?;

    let params: Params = match Params::from_json_sting(&fs::read_to_string(param)?) {
        Ok(params) => params,
        Err(e) => {
            log_line(logger, &format!("Error reading param file: {}", e))?;
            tcs_report.add_error(e.to_string());
            return Ok((tcs_report, Some((r1_file.clone(), r2_file.clone()))));
        }
    };

    // Validate the params
    // This will validate the params and return a validated Params struct.
    // If the params are not valid, it will return an error.
    // The validated Params struct will be used to filter the R1 and R2 pairs.
    // If there is an error validating the params, it will log the error to the logger and return a TcsReport with the error.
    // The TcsReport with error will be handled in the downstream processing.
    log_line(logger, "Validating Params")?;

    let validated_params = match params.validate() {
        Ok(validated_params) => validated_params,
        Err(e) => {
            log_line(logger, &format!("Error validating params: {}", e))?;
            tcs_report.add_error(e.to_string());
            return Ok((tcs_report, Some((r1_file.clone(), r2_file.clone()))));
        }
    };

    let pairs = match read_fastq_file(&fastq_files) {
        Ok(pairs) => pairs,
        Err(e) => {
            log_line(logger, &format!("Error reading fastq files: {}", e))?;
            tcs_report.add_error(e.to_string());
            return Ok((tcs_report, Some((r1_file.clone(), r2_file.clone()))));
        }
    };

    log_line(logger, "Reading Fastq files")?;
    log_line(
        logger,
        &format!("Number of raw fastq records: {}", pairs.len()),
    )?;

    tcs_report.set_total_reads(pairs.len());

    // Process the pairs in parallel
    // This will filter the R1 and R2 pairs based on the validated params.
    // It will use Rayon to process the pairs in parallel.
    // The filter_r1_r2_pairs function will return a PairedRecordFilterResult enum.
    // If the pair is valid, it will return a FilteredPair struct.
    // If the pair is invalid, it will return a reason for failure.
    // If there is an error processing the pairs, it will log the error to the logger and return a TcsReport with the error.
    // The TcsReport with error will be handled in the downstream processing.
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

    // log the de-multiplexed pairs
    // these de-muliplexed pairs are grouped by region, and will be processed downstream.
    // we log the number of valid pairs for each region.
    log_line(logger, "De-multiplexed pairs")?;
    for (region, filtered_pairs) in &groups {
        log_line(
            logger,
            &format!(
                "Region: {}, valid r1 r2 pairs: {}",
                region,
                filtered_pairs.len()
            ),
        )?;
    }

    // for the failed pairs, we log the total number of failed pairs and the reasons for failure.
    // but we do not log the individual pairs in the log file. We will populate a summary of the reasons for failure as part of the output.
    log_line(logger, "Failed pairs")?;
    let mut fail_frequency = HashMap::new();
    for fail in &fails {
        *fail_frequency.entry(fail.to_string()).or_insert(0) += 1;
        tcs_report.add_failed_match_reason(fail.to_string());
    }

    log_line(
        logger,
        &format!(
            "A total of {} paired sequences failed to map to de-multiplex for {} number of reasons",
            fails.len(),
            fail_frequency.len()
        ),
    )?;

    // Log all the reasons of errors. Errors are not expected to be many, so we log them all.
    // Also errors are different from fails, errors are unexpected issues that occur during processing.

    if errors.is_empty() {
        log_line(logger, "No errors encountered when filtering raw sequences")?;
    } else {
        log_line(
            logger,
            &format!(
                "A total of {} errors encountered when filtering raw sequences",
                errors.len()
            ),
        )?;
    }

    // log each error in the tcs_report and the logger.
    // These errors will not stop the processing, but will be logged for debugging purposes in the Warning section of the report.
    // We will also add these errors to the TcsReportWarnings enum, with a type of R1R2filteringwarning.
    for error in errors {
        log_line(logger, &format!("Error: {}", error))?;
        tcs_report.add_warning(TcsReportWarnings::R1R2filteringwarning(error.to_string()));
    }

    // now processing consensus calling
    // We will use the ConsensusStrategy::Weighted with the steepness and midpoint parameters.
    // The ConsensusStrategy::Weighted will use the steepness and midpoint parameters to calculate the consensus sequence.
    // The steepness parameter will control the steepness of the curve, and the midpoint parameter will control the midpoint of the curve.
    // The ConsensusStrategy::Weighted will be used to calculate the consensus sequence for each region
    let consensus_strategy =
        ConsensusStrategy::Weighted(ConsensusParams::new(steepness as f64, midpoint as f64));

    log_line(logger, "Starting consensus calling")?;

    // Process each region in sequence
    // We will iterate over each region and call the TcsConsensus::build_from_filtered_pairs function to build the consensus sequence for each region.
    // The TcsConsensus::build_from_filtered_pairs function will take the filtered pairs for each region and return a tuple containing the consensus results, errors, and UMI summary.
    // If there is an error processing a region, we will log the error and continue to the next region.
    // We will also create a RegionReport for each region and add it to the TcsReport.
    let mut region_reports = Vec::new(); // This will hold the reports for each region for the field `region_reports` in TcsReport
    for (region, filtered_pairs) in &groups {
        let region_params =
            validated_params
                .get_region_params(region)
                .ok_or(TcsError::UnexpectedError(format!(
                    "No parameters found for region: {}",
                    region
                )))?;
        let mut region_report = RegionReport::new();
        region_report.set_region_name(region.clone());
        region_report.set_filtered_reads_for_region(filtered_pairs.len());
        log_line(
            logger,
            &format!(
                "Processing region: {}, with {} valid pairs",
                region,
                filtered_pairs.len()
            ),
        )?;

        // Build consensus for the region
        // This will call the TcsConsensus::build_from_filtered_pairs function to build the consensus sequence for the region.
        // We use Rayon to process the pairs in parallel.
        // We capture UMIDistError, if it occurs, we log it in the logger and add a warning to the TcsReport, and continue to the next region.
        // If the consensus calling is successful, we will log the number of UMIs found, the UMI cut-off, and the number of UMIs passing the error cutoff.
        // TcsConsensus will be collected into a vector of TcsConsensus, which will be stored in the region_report.
        // Errors during consensus calling for individual UMI families will be logged, and warnings will be added to the TcsReport.
        // The UMI summary will be collected and added to the RegionReport as part of the TcsReport.

        let (mut consensus_results, consensus_errors, umi_summary) =
            match TcsConsensus::build_from_filtered_pairs(
                filtered_pairs,
                consensus_strategy,
                params.platform_error_rate,
            ) {
                Ok(tcs_consensus_building_output) => (
                    tcs_consensus_building_output.tcs_consensus().clone(),
                    tcs_consensus_building_output.errors().clone(),
                    tcs_consensus_building_output.umi_summary().clone(),
                ),
                Err(e) => {
                    log_line(
                        logger,
                        &format!("UMI Distribution Error for Region {}: {}", region, e),
                    )?;
                    tcs_report.add_warning(TcsReportWarnings::UMIDistErrorWithRegion(
                        region.clone(),
                        e.to_string(),
                    ));
                    region_reports.push(region_report);
                    continue; // Skip to the next region if there's an error
                }
            };

        for err in consensus_errors {
            tcs_report.add_warning(TcsReportWarnings::ConsensusErrorIndividualWithRegion(
                region.clone(),
                err.clone(),
            ));

            log_line(
                logger,
                &format!("Consensus Error for Region {}: {}", region, err),
            )?;
        }

        let passed_umi_families_distribution = umi_summary.get_passed_umis_hashmap();
        log_line(
            logger,
            &format!(
                "Region: {}, A total of {} UMIs found, UMI cut-off is {}, a total of {} UMIs passing the error cutoff",
                region,
                umi_summary.umi_freq().len(),
                umi_summary.umi_cut_off(),
                passed_umi_families_distribution.len()
            ),
        )?;
        log_line(
            logger,
            &format!(
                "Number of consensus sequences generated for region {}: {}",
                region,
                consensus_results.len()
            ),
        )?;

        // Start end-joining for the region

        if region_params.end_join {
            log_line(logger, &format!("End-joining for region: {}", region))?;
        } else {
            log_line(
                logger,
                &format!(
                    "End-joining not required for {}, skip end-joining, QC and Trimming",
                    region
                ),
            )?;
            region_report.set_tcs_consensus_results(Some(consensus_results));
            region_reports.push(region_report);
            continue; // Continue to the next region if end-joining is not required
        }

        // End-joining logic below
        if consensus_results.is_empty() {
            log_line(
                logger,
                &format!(
                    "No consensus sequences generated for region: {}. Skipping end-joining.",
                    region
                ),
            )?;
            region_report.set_tcs_consensus_results(Some(consensus_results));
            region_reports.push(region_report);
            continue; // Skip end-joining if no consensus sequences are available
        }
        if let Err(error) = join_consensus_fastq_vec(
            &mut consensus_results,
            region_params.end_join_option,
            region_params.overlap as usize,
        ) {
            log_line(
                logger,
                &format!(
                    "Error during end-joining for region: {}. Invidual consensus sequences will not be end-joined. Error: {}",
                    region, error
                ),
            )?;
            tcs_report.add_warning(TcsReportWarnings::EndJoiningErrorWithRegion(
                region.clone(),
                error.to_string(),
            ));
        };

        log_line(
            logger,
            &format!(
                "End-joining completed for region: {} without warnings",
                region
            ),
        )?;

        // TODO: QC and trimming logic

        if region_params.tcs_qc {
            log_line(logger, &format!("QC (and trimming) for region: {}", region))?;

            if let Err(error) = qc_and_trim_consensus_fastq_vec(
                &mut consensus_results,
                region_params.qc_config.as_ref(),
                region_params.trim_config.as_ref(),
            ) {
                log_line(
                    logger,
                    &format!(
                        "Error during QC and trimming for region: {}. Error: {}",
                        region, error
                    ),
                )?;
                tcs_report.add_warning(TcsReportWarnings::QcAndTrimErrorWithRegion(
                    region.clone(),
                    error.to_string(),
                ));
            } else {
                log_line(
                    logger,
                    &format!(
                        "QC and trimming completed for region: {}, a total of {} QC/Trimmed TCS obtained",
                        region,
                        count_passed(&consensus_results)
                    ),
                )?;
            }
        } else {
            log_line(
                logger,
                &format!("QC not required for {}, skip QC and Trimming", region),
            )?;
            for consensus_result in &mut consensus_results {
                // If QC is not required, we still need to set the trimmed flag to false
                consensus_result.set_qc(TcsConsensusQcResult::NotRequired);
            }
            region_report.set_tcs_consensus_results(Some(consensus_results));
            region_reports.push(region_report);
        }
    }

    // TODO: downstream processing
    // 1. consensus calling for each region. DONE!
    // 2. end-joining of consensus sequences. DONE!
    // 3. qc
    // 4. trimming
    // 5. write fastq and fasta files.
    // 6. write UMI files in JSON format.
    // 7. export a summary report.
    // 8. Error handling and logging. Some errors are expected, so we do not panic, but log them and continue processing.

    Ok((tcs_report, Some((r1_file.clone(), r2_file.clone()))))
}

//TODO: write details of the function
pub fn tcs_write(
    tcs_report: &TcsReport,
    logger: &mut BufWriter<File>,
) -> Result<(), Box<dyn Error>> {
    todo!(
        "TCS write function called, but not implemented yet. This function should write the TCS report to a file. with the following details: {:?} and logger: {:?}",
        tcs_report,
        logger
    )
}

fn tcs_init(input: &str) -> Result<(TcsReport, BufWriter<File>), Box<dyn Error>> {
    // Initialize the TCS report
    let mut tcs_report = TcsReport::new();

    let input_dir = Path::new(input);
    // Check if the input directory exists
    if !input_dir.exists() {
        return Err(TcsError::InputDirNotFound(input.to_string()).into());
    }
    // Check if the input directory is a valid directory
    if !input_dir.is_dir() {
        return Err(TcsError::NotADirectory(input.to_string()).into());
    }

    tcs_report.set_input_directory(input.to_string());

    let logfile = OpenOptions::new()
        .create(true)
        .append(true)
        .open(input_dir.join("run_log.txt"))?;

    let logger: BufWriter<File> = BufWriter::new(logfile);

    Ok((tcs_report, logger))
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
