use std::fs::File;
use std::io::BufReader;

use bio::io::fastq::{self, Record};
use flate2::read::MultiGzDecoder;

use crate::utils::tcs_helper::{DataType, FastqFiles};

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
