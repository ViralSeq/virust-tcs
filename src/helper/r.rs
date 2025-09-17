use std::error::Error;
use std::process::Command;

pub fn check_r_installed() -> Result<(), Box<dyn Error>> {
    // In a real implementation, you would check if R is available in the system PATH
    let output = Command::new("R").arg("--version").output();

    match output {
        Ok(out) => {
            if out.status.success() {
                Ok(())
            } else {
                Err("R command failed to execute properly.".into())
            }
        }
        Err(_) => Err("R is not installed or not found in PATH.".into()),
    }
}

pub fn check_r_packages() -> Result<(), Box<dyn Error>> {
    let check_env_script = include_str!("../../resources/r_scripts/check_env.r");

    let output = Command::new("Rscript")
        .arg("-e")
        .arg(check_env_script)
        .output()?;

    if output.status.success() {
        Ok(())
    } else {
        let stderr = String::from_utf8_lossy(&output.stderr);
        Err(format!("R packages check failed: {}", stderr).into())
    }
}

pub fn get_sdrm_r_script() -> &'static str {
    include_str!("../../resources/r_scripts/sdrm_r.r")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_check_r_installed() {
        // This test will pass if R is installed on the system where the test is run.
        // In a real unit test, you might want to mock the Command::new call.
        match check_r_installed() {
            Ok(_) => println!("Found R installation."),
            Err(e) => println!("R is not installed: {}", e),
        }
    }

    #[test]
    fn test_check_r_packages() {
        match check_r_packages() {
            Ok(_) => println!("R packages are installed."),
            Err(e) => println!("R packages check failed: {}", e),
        }
    }

    #[test]
    fn test_get_sdrm_r_script() {
        let script = get_sdrm_r_script();
        assert_eq!(script.len(), 1080); // Adjust this length if the script changes
    }
}
