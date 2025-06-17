pub mod error;
pub mod fastq_files;
pub mod filter_r1_r2;
pub mod tcs_consensus;
pub mod tcs_qc;
pub mod tcs_report;
pub mod tcs_summary;
pub mod utils;

pub use error::TcsError;
pub use fastq_files::{FastqFiles, validate_files};
pub use filter_r1_r2::{FilteredPair, PairedRecordFilterResult, filter_r1_r2_pairs};
pub use tcs_consensus::{
    TcsConsensus, TcsConsensusBuildingOutput, TcsConsensusQcResult, count_passed,
    join_consensus_fastq_vec, qc_and_trim_consensus_fastq_vec,
};
pub use tcs_qc::{QcAlgorithm, QcReference, TcsQcInput};
pub use tcs_report::{AdvancedSettings, RegionReport, TcsReport, TcsReportWarnings};
pub use utils::*;
