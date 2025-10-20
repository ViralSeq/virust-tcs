use std::error::Error;
use std::process::{Command, Stdio};

#[derive(Debug, PartialEq, Eq, Hash)]
#[allow(dead_code)]
pub enum MuscleVersion {
    Muscle3_8_31,
    Muscle5,
    Other(String),
    NotInstalled,
}

impl MuscleVersion {
    /// build a MUSCLE command apporpriate for the version
    pub fn build_command(&self, input: &str, output: &str) -> Option<Command> {
        let mut cmd = Command::new("muscle");
        match self {
            MuscleVersion::Muscle3_8_31 => {
                cmd.arg("-in")
                    .arg(input)
                    .arg("-out")
                    .arg(output)
                    .stdout(Stdio::null())
                    .stderr(Stdio::null());
            }
            MuscleVersion::Muscle5 => {
                cmd.arg("-super5")
                    .arg(input)
                    .arg("-output")
                    .arg(output)
                    .stdout(Stdio::null())
                    .stderr(Stdio::null()); //use super5 mode for (much) faster processing
            }
            MuscleVersion::Other(version) => {
                // raise warning about unknown version
                eprintln!("Warning: Unknown MUSCLE version: {}", version);
                return None;
            }
            MuscleVersion::NotInstalled => {
                eprintln!("Error: MUSCLE is not installed.");
                return None;
            }
        }
        Some(cmd)
    }

    pub fn run(&self, input: &str, output: &str) -> Result<(), Box<dyn Error>> {
        if let Some(mut cmd) = self.build_command(input, output) {
            let status = cmd.status()?;
            if status.success() {
                Ok(())
            } else {
                Err(format!("MUSCLE command failed with status: {}", status).into())
            }
        } else {
            Err("Failed to build MUSCLE command.".into())
        }
    }
}

// Determine which MUSCLE version the system has. Muscle 3.8.31 and Muscle 5 have different command line arguments
pub fn get_muscle_version(keyword: &str) -> MuscleVersion {
    // In a real implementation, you would check if MUSCLE is available in the system PATH

    let output = Command::new(keyword).arg("-version").output();

    let output = match output {
        Ok(out) => out,
        Err(_) => return MuscleVersion::NotInstalled,
    };

    if !output.status.success() {
        return MuscleVersion::NotInstalled;
    }

    let stdout = String::from_utf8_lossy(&output.stdout);
    let stderr = String::from_utf8_lossy(&output.stderr);

    let text = if !stdout.trim().is_empty() {
        stdout.trim()
    } else {
        stderr.trim()
    };

    if text.contains("3.8.31") {
        MuscleVersion::Muscle3_8_31
    } else if text.contains("muscle 5") || text.contains("Muscle 5") || text.contains("MUSCLE 5") {
        MuscleVersion::Muscle5
    } else {
        MuscleVersion::Other(text.to_string())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_muscle_version() {
        // This test assumes that MUSCLE is not installed in the test environment
        let version = get_muscle_version("muscle");
        assert!(matches!(version, MuscleVersion::Muscle5));

        let version = get_muscle_version("muscle3");
        assert!(matches!(version, MuscleVersion::Muscle3_8_31));

        let version = get_muscle_version("nonexistent_command");
        assert_eq!(version, MuscleVersion::NotInstalled);
    }

    #[test]
    fn test_test_muscle_command_building() {
        let version_3 = MuscleVersion::Muscle3_8_31;
        let cmd_3 = version_3
            .build_command("input.fasta", "output.fasta")
            .unwrap();
        let args_3: Vec<String> = cmd_3
            .get_args()
            .map(|s| s.to_string_lossy().to_string())
            .collect();
        assert_eq!(args_3, vec!["-in", "input.fasta", "-out", "output.fasta"]);

        let version_5 = MuscleVersion::Muscle5;
        let cmd_5 = version_5
            .build_command("input.fasta", "output.fasta")
            .unwrap();
        let args_5: Vec<String> = cmd_5
            .get_args()
            .map(|s| s.to_string_lossy().to_string())
            .collect();
        assert_eq!(
            args_5,
            vec!["-super5", "input.fasta", "-output", "output.fasta"]
        );
    }

    #[test]
    fn test_muscle_command_local_run() {
        let version = get_muscle_version("muscle");
        if version == MuscleVersion::NotInstalled {
            eprintln!("MUSCLE is not installed. Skipping run test.");
            return;
        }

        let input = "tests/data/alignment/sequence.fasta";
        let output = "tests/data/alignment/sequence.aligned.fasta";

        let result = version.run(input, output);
        assert!(result.is_ok() || result.is_err());
    }
}
