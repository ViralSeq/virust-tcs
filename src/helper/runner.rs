#![allow(dead_code)]

use std::error::Error;
use std::process::Command;

const QMD_INDEX_TEMPLATE: &str = "resources/quarto/index.qmd";
const QMD_SAMPLE_TEMPLATE: &str = "resources/quarto/sample_report.qmd";
const R_CHECK_ENV: &str = "resources/r_scripts/check_env.R";

/// Check if R, Quarto, and Python3 are installed, necessary for report generation.
pub fn check_r_installed() -> Result<(), Box<dyn Error>> {
    is_available("Rscript", &["--version"], "R")
}

pub fn check_quarto_installed() -> Result<(), Box<dyn Error>> {
    is_available("quarto", &["--version"], "Quarto")
}

pub fn check_python3_installed() -> Result<(), Box<dyn Error>> {
    is_available("python3", &["--version"], "Python 3")
}

fn is_available(cmd: &str, args: &[&str], program: &str) -> Result<(), Box<dyn Error>> {
    Command::new(cmd)
        .args(args)
        .output()
        .map(|output| {
            if output.status.success() {
                Ok(())
            } else {
                Err(format!("{} command failed to execute properly.", cmd).into())
            }
        })
        .unwrap_or_else(|_| {
            Err(format!("{} is not installed or not found in PATH.", program).into())
        })
}
