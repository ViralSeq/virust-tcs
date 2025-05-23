pub mod consensus;
pub mod io;
pub mod params;
pub mod tcs_helper;
pub mod umi;
pub mod umis;

use thiserror::Error;

//TODO: write details of the enum
#[derive(Error, Debug)]
pub enum ParamValidationError {
    #[error("cDNA primer is not valid IUPAC sequence: {0}")]
    CDNAprimerInvalid(String),

    #[error("Forward primer is not valid IUPAC sequence: {0}")]
    ForwardPrimerInvalid(String),

    #[error("UMI parsing failed: {0}")]
    UmiError(#[from] UmiError),
}

#[derive(Error, Debug)]
pub enum UmiError {
    #[error("More than one regular UMI found in cDNA primer: {0}")]
    DuplicatedRegularUMI(String),

    #[error("No UMI found in cDNA primer: {0}")]
    NoUMIFound(String),
}
