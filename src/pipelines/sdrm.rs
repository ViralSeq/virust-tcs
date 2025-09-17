//TODO: SDRM pipeline

use std::error::Error;

use crate::helper::muscle::get_muscle_version;
use crate::helper::r::{check_r_installed, get_sdrm_r_script};

pub fn run_sdrm(input: String, version: String) -> Result<(), Box<dyn Error>> {
    // Placeholder implementation
    println!(
        "Running SDRM pipeline with input: {}, version: {}",
        input, version
    );

    // check environment, ensure MSA aligner (MUSCLE) is available

    let muscle_version = get_muscle_version("muscle");

    println!("Detected MUSCLE version: {:?}", muscle_version); //placeholder

    // check if R and required R packages are installed

    check_r_installed()?;

    let r_script: &'static str = get_sdrm_r_script();

    println!("Using R script:\n{}", r_script); //placeholder

    todo!("Implement the SDRM pipeline logic here");
}
