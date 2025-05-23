use chrono::Local;
use regex::Regex;
use std::error::Error;
use std::fs;
use std::fs::File;
use std::io::BufWriter;
use std::io::{Result as IoResult, Write};
use std::path::{Path, PathBuf};

#[derive(Debug, PartialEq)]
pub enum DataType {
    Fastq,
    FastqGz,
}

#[derive(Debug)]
pub struct FastqFiles {
    pub r1_file: PathBuf,
    pub r2_file: PathBuf,
    pub data_type: DataType,
}

pub fn validate_files(input: &str) -> Result<FastqFiles, Box<dyn Error>> {
    // Check if the input file exists
    if !Path::new(input).exists() {
        return Err(format!("Input directory {} does not exist", input).into());
    }

    // Check if the input file is a valid file
    if !Path::new(input).is_dir() {
        return Err(format!("Input path {} is not a valid directory", input).into());
    }

    let entries = fs::read_dir(input)?;
    let mut r1_candidates = vec![];
    let mut r2_candidates = vec![];

    let r1_re = Regex::new(r"(?i)(^|[_\-.])r1([_\-.]\d+)?\.f(ast)?q(\.gz)?$")?;
    let r2_re = Regex::new(r"(?i)(^|[_\-.])r2([_\-.]\d+)?\.f(ast)?q(\.gz)?$")?;

    for entry in entries {
        let entry = entry?;
        let path = entry.path();
        if let Some(fname) = path.file_name().and_then(|f| f.to_str()) {
            if r1_re.is_match(fname) {
                r1_candidates.push(path);
            } else if r2_re.is_match(fname) {
                r2_candidates.push(path);
            }
        }
    }

    // Error: check number of files
    match (r1_candidates.len(), r2_candidates.len()) {
        (0, 0) => {
            return Err("No R1 or R2 files found in the input directory".into());
        }
        (0, _) => {
            return Err("No R1 files found in the input directory".into());
        }
        (_, 0) => {
            return Err("No R2 files found in the input directory".into());
        }
        (1, 1) => {
            // Do nothing, valid case
        }
        (n, m) => {
            return Err(format!(
                "Found {} R1 files and {} R2 files. Expected 1 of each.",
                n, m
            )
            .into());
        }
    }

    let r1_file = &r1_candidates[0];
    let r2_file = &r2_candidates[0];

    let r1_gz = r1_file.extension().map(|ext| ext == "gz").unwrap_or(false);
    let r2_gz = r2_file.extension().map(|ext| ext == "gz").unwrap_or(false);

    // check type consistency
    if r1_gz != r2_gz {
        return Err(format!(
            "File type mismatch: R1 is {}compressed, R2 is {}compressed.",
            if r1_gz { "" } else { "not " },
            if r2_gz { "" } else { "not " },
        )
        .into());
    }

    let data_type = if r1_gz {
        DataType::FastqGz
    } else {
        DataType::Fastq
    };
    let fastq_files = FastqFiles {
        r1_file: r1_file.clone(),
        r2_file: r2_file.clone(),
        data_type,
    };

    Ok(fastq_files)
}

pub fn log_line(writer: &mut BufWriter<File>, message: &str) -> IoResult<()> {
    let now = Local::now().format("%Y-%m-%d %H:%M:%S");
    writeln!(writer, "[{}] {}", now, message)?;
    writer.flush()?;
    Ok(())
}

pub fn diff_positions(a: &str, b: &str) -> Vec<usize> {
    a.chars()
        .zip(b.chars())
        .enumerate()
        .filter_map(|(i, (ca, cb))| if ca != cb { Some(i) } else { None })
        .collect()
}

pub fn diff_byte_equal_length(a: &[u8], b: &[u8]) -> Vec<usize> {
    (0..a.len()).filter(|&i| a[i] != b[i]).collect()
}

#[cfg(test)]

mod tests {
    use super::*;
    use std::path::PathBuf;

    #[test]
    fn test_validate_files() {
        let input = "tests/data/hivdr_control";
        let result = validate_files(input);
        assert!(result.is_ok());
        let fastq_files = result.unwrap();
        assert_eq!(
            fastq_files.r1_file,
            PathBuf::from("tests/data/hivdr_control/r1.fastq.gz")
        );
        assert_eq!(
            fastq_files.r2_file,
            PathBuf::from("tests/data/hivdr_control/r2.fastq.gz")
        );
        assert_eq!(fastq_files.data_type, DataType::FastqGz);

        let input = "tests/data/some_non_existent_directory";
        let result = validate_files(input);
        assert!(result.is_err());

        let input = "tests/data/test_dir";
        let result = validate_files(input);
        assert!(result.is_ok());
        let fastq_files = result.unwrap();
        assert_eq!(
            fastq_files.r1_file,
            PathBuf::from("tests/data/test_dir/mydata_R1_001.fastq")
        );
        assert_eq!(
            fastq_files.r2_file,
            PathBuf::from("tests/data/test_dir/mydata_R2_001.fastq")
        );
        assert_eq!(fastq_files.data_type, DataType::Fastq);

        let input = "tests/data/test_dir2";
        let result = validate_files(input);
        assert!(result.is_err());
        assert!(
            result.unwrap_err().to_string()
                == "Found 2 R1 files and 1 R2 files. Expected 1 of each."
        );
    }
}
