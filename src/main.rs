use clap::Parser;
use virust_tcs::cli::Args;
use virust_tcs::cli::Commands;
use virust_tcs::params_generator;

fn main() {
    let args = Args::parse();

    match args.command {
        Commands::Run {
            input,
            param,
            keep_original,
        } => {
            println!(
                "Running TCS pipeline with input: {}, param: {}, keep_original: {}",
                input, param, keep_original
            );
            // TODO: Call the function to run the pipeline here
            todo!();
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
            println!(
                "Running TCS DR pipeline with input: {}, version: {}, keep_original: {}",
                input, version, keep_original
            );
            // TODO: Call the function to run the DR pipeline here
            todo!();
        }
        Commands::DrParams { version } => {
            println!("Listing DR params...");
            if let Some(v) = version {
                println!("Version: {}", v);
            } else {
                println!("Listing all available versions...");
            }
            // TODO: Call the function to list DR params here
            todo!()
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
