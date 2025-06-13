use serde::{Deserialize, Serialize};
use thiserror::Error;

#[derive(Error, Debug, Clone, Serialize, Deserialize)]
pub enum TcsError {
    #[error("Input directory does not exist: {0}")]
    InputDirNotFound(String),
    #[error("Input path is not a valid directory: {0}")]
    NotADirectory(String),
    #[error("No R1 or R2 files found in the input directory")]
    NoFastqFilesFound,
    #[error("No R1 files found in the input directory")]
    NoR1FilesFound,
    #[error("No R2 files found in the input directory")]
    NoR2FilesFound,
    #[error("Found {0} R1 files and {1} R2 files. Expected 1 of each.")]
    MultipleFilesFound(usize, usize),
    #[error("File type mismatch: R1 is {0}compressed, R2 is {1}compressed.")]
    FileTypeMismatch(String, String),
    #[error("Invalid R1 header: {0}")]
    InvalidR1Header(String),
    #[error("Invalid R2 header: {0}")]
    InvalidR2Header(String),
    #[error("Empty fastq record")]
    EmptyFastqRecord,
    #[error("R1 R2 header mismatch: R1: {0}, R2: {1}")]
    R1R2HeaderMismatch(String, String),
    #[error("Invalid R1 record: {0}")]
    InvalidR1Record(String),
    #[error("Invalid R2 record: {0}")]
    InvalidR2Record(String),
    #[error(
        "Invalid read length: Platform Format: {0}, should be equal or less to Read 1 Length: {1} and Read 2: {2}"
    )]
    InvalidReadLength(usize, usize, usize),
    #[error("Failed to access the param file from the given path: {0}")]
    ParamFileAccessError(String),
    #[error("Unexpected error: {0}")]
    UnexpectedError(String),
}
