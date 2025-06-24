use clap::Parser;
use virust_tcs::cli::Args;
use virust_tcs::cli::Commands;
use virust_tcs::helper::*;
use virust_tcs::pipelines::params_generator;
use virust_tcs::pipelines::tcs::*;

fn main() {
    let args = Args::parse();

    match args.command {
        Commands::Run {
            input,
            param,
            keep_original,
            steepness,
            midpoint,
        } => {
            let params_input_type = ParamsInputType::FromFilePath(param.clone());
            tcs(
                &input,
                params_input_type,
                keep_original,
                steepness,
                midpoint,
            )
            .unwrap_or_else(|err| {
                eprintln!("Fatal Error: {} occurred during processing", err);
                std::process::exit(1);
            });
        }
        Commands::Generate {} => {
            // Call the function to generate the param file here
            params_generator::exec();
        }
        Commands::DR {
            input,
            version,
            keep_original,
        } => {
            let params_input_type = ParamsInputType::PresetID(version.clone().to_lowercase());
            tcs(
                &input,
                params_input_type,
                keep_original,
                consensus::DEFAULT_K as f32,
                consensus::DEFAULT_Q0 as u8,
            )
            .unwrap_or_else(|err| {
                eprintln!("Fatal Error: {} occurred during processing", err);
                std::process::exit(1);
            });
        }
        Commands::DrParams { version } => {
            println!("Listing DR params...");
            let all_versions = params::dr_presets_names();
            if let Some(v) = version {
                let params_query = params::Params::from_preset(&v.to_lowercase());
                if params_query.is_err() {
                    eprintln!("Error: {}", params_query.unwrap_err());
                    std::process::exit(1);
                } else {
                    println!("Params for version {}:\n{}", v, params_query.unwrap());
                }
            } else {
                println!(
                    "Available versions for DR params: {}",
                    all_versions.join(", ")
                );
            }
        }
        Commands::SDRM { input, version } => {
            println!(
                "Running SDRM pipeline with input: {}, version: {}",
                input, version
            );
            // TODO: Call the function to run the SDRM pipeline here
            todo!();
        }
        Commands::Log { input } => {
            println!("Running TCS log pipeline with input: {}", input);
            // TODO: Call the function to run the log pipeline here
            todo!();
        }
    }
}
