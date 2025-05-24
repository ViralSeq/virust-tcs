use std::fs::File;
use std::io::BufReader;

use bio::io::fastq::{self, Record};
use flate2::read::MultiGzDecoder;

use crate::utils::tcs_helper::{DataType, FastqFiles};

/// Reads paried R1 R2 fastq files and returns a vector of tuples containing the records from both files.
/// The function takes a `FastqFiles` struct as an argument, which contains the paths to the R1 and R2 files.
/// The function uses the `bio` crate to read the fastq files and returns a vector of tuples containing the records from both files.
/// The function also handles different data types (Fastq and FastqGz) using the `DataType` enum.
/// The function returns a `Result` containing a vector of tuples of records or an `io::Error` if there was an error reading the files.
/// The function uses flate2 to handle gzipped files.
/// *Arguments*
/// - `files`: A `FastqFiles` struct containing the paths to the R1 and R2 files.
/// *Returns*
/// - `Result<Vec<(Record, Record)>, std::io::Error>`: A result containing a vector of tuples of records or an `io::Error` if there was an error reading the files.
pub fn read_fastq_file(files: &FastqFiles) -> std::io::Result<Vec<(Record, Record)>> {
    let r1_file = File::open(&files.r1_file)?;
    let r2_file = File::open(&files.r2_file)?;

    let (r1_stream, r2_stream): (Box<dyn std::io::Read>, Box<dyn std::io::Read>) =
        match files.data_type {
            DataType::Fastq => (
                Box::new(BufReader::new(r1_file)),
                Box::new(BufReader::new(r2_file)),
            ),
            DataType::FastqGz => (
                Box::new(MultiGzDecoder::new(BufReader::new(r1_file))),
                Box::new(MultiGzDecoder::new(BufReader::new(r2_file))),
            ),
        };
    let r1_reader = fastq::Reader::new(BufReader::new(r1_stream));
    let r2_reader = fastq::Reader::new(BufReader::new(r2_stream));

    // Collect record pairs into Vec
    let pairs: Vec<(Record, Record)> = r1_reader
        .records()
        .zip(r2_reader.records())
        .filter_map(|(r1, r2)| match (r1.ok(), r2.ok()) {
            (Some(rec1), Some(rec2)) => Some((rec1, rec2)),
            _ => None,
        })
        .collect();
    Ok(pairs)
}
