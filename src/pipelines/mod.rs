use crate::utils::ParamValidationError;
use thiserror::Error;

pub mod log;
pub mod params_generator;
pub mod sdrm;
pub mod tcs;

//TODO:  write details of the enum
#[derive(Error, Debug)]
pub enum PipelineError {}

//TODO: continue writing details of the enum
#[derive(Error, Debug)]
pub enum TCSError {
    #[error("Param Validation Error: {0}")]
    ParamValidationError(#[from] ParamValidationError),
}
