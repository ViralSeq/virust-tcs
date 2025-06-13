use crate::cli::BANNER;
use crate::helper::params::Params;
use crate::helper::params::RegionParams;
use crate::helper::params::{validate_cdna_primer, validate_nt_words};
use std::fs::File;
use std::io::{self, Write};
use std::path::PathBuf;

pub fn exec() {
    println!("{}", BANNER);

    println!("{}", "-".repeat(58));
    println!(
        "| JSON Parameter Generator for TCS {} by Shuntai Zhou |",
        env!("CARGO_PKG_VERSION")
    );
    println!("{}", "-".repeat(58));

    print!(
        "Enter the path to the directory that contains the MiSeq pair-end R1 and R2 .fastq or .fastq.gz files (e.g. /path/to/fastq_dir):\n>  "
    );

    let input_dir = collect_input();

    print!("Choose MiSeq Platform (1-3):\n1. 150x7x150\n2. 250x7x250\n3. 300x7x300 (default)\n>  ");

    let platform = match collect_input().as_str() {
        "1" => 150,
        "2" => 250,
        _ => 300,
    };

    print!(
        "Enter the estimated platform error rate (for TCS cut-off calculation), default as 0.02:\n>  "
    );
    let error_rate = match collect_input().as_str() {
        "" => 0.02,
        input => input.parse::<f32>().unwrap_or(0.02),
    };

    let mut regions: Vec<RegionParams> = Vec::new();

    loop {
        let region_name = loop {
            print!("Enter the name for the sequenced region: \n>  ");
            let region_name = collect_input();
            if region_name.is_empty() {
                eprintln!("Region name cannot be empty. Please enter a valid name.");
                continue;
            }
            if region_name.len() > 20 {
                eprintln!(
                    "Region name is too long. Please enter a name with 20 characters or less.\n"
                );
                continue;
            } else {
                break region_name;
            }
        };

        let cdna_primer = loop {
            print!("Enter the cDNA primer sequence:\n>  ");
            let cdna_primer = collect_input().to_uppercase();

            if validate_cdna_primer(&cdna_primer).is_ok() {
                break cdna_primer;
            } else {
                eprintln!("Invalid cDNA primer sequence. Please enter a valid sequence.");
                continue;
            }
        };

        let forward_primer = loop {
            print!("Enter the forward primer sequence:\n>  ");
            let forward_primer = collect_input().to_uppercase();

            if validate_nt_words(&forward_primer).is_ok() {
                break forward_primer;
            } else {
                eprintln!("Invalid forward primer sequence. Please enter a valid sequence.");
                continue;
            }
        };

        print!("Enter supermajority cut-off (0.5 - 1.0). Default Simple Majority (0.5):\n>  ");
        let mut majority_cutoff = match collect_input().as_str() {
            "" => 0.5,
            input => input.parse::<f32>().unwrap_or(0.5),
        };

        majority_cutoff = if majority_cutoff < 0.5 {
            0.5
        } else if majority_cutoff > 1.0 {
            1.0
        } else {
            majority_cutoff
        };

        print!("Need end-joining? (y/n, default as n):\n>  ");

        let end_join = match collect_input().as_str() {
            "y" | "Y" => true,
            _ => false,
        };

        let end_join_option: u32 = if end_join {
            collect_end_join_option() as u32
        } else {
            0
        };

        let overlap_size: u32 = if end_join_option == 2 {
            print!("Enter the overlap size (default as 20):\n>  ");
            match collect_input().as_str() {
                "" => 20,
                input => input.parse::<u32>().unwrap_or(20),
            }
        } else {
            0
        };

        print!("Need TCS QC? Support HIV-1 and SVI (y/n, default as n):\n>  ");
        let tcs_qc = match collect_input().as_str() {
            "y" | "Y" => true,
            _ => false,
        };

        let (ref_genome, ref_start, ref_start_lower, ref_end, ref_end_lower) = if tcs_qc {
            get_ref_and_locations()
        } else {
            (String::new(), 0, None, 0, None)
        };

        let indel = if tcs_qc {
            print!("allow indels in QC (y/n, default as n):\n>  ");
            match collect_input().as_str() {
                "y" | "Y" => true,
                _ => false,
            }
        } else {
            false
        };

        print!("Need trimming? (y/n, default as n):\n>  ");
        let trim = match collect_input().as_str() {
            "y" | "Y" => true,
            _ => false,
        };

        let (trim_ref, trim_ref_start, _, trim_ref_end, _) = if trim {
            get_ref_and_locations()
        } else {
            (String::new(), 0, None, 0, None)
        };
        let trim_ref = if trim_ref.is_empty() {
            None
        } else {
            Some(trim_ref)
        };
        let trim_ref_start = if trim_ref_start == 0 {
            None
        } else {
            Some(trim_ref_start)
        };
        let trim_ref_end = if trim_ref_end == 0 {
            None
        } else {
            Some(trim_ref_end)
        };
        regions.push(RegionParams {
            region: region_name,
            forward: forward_primer,
            cdna: cdna_primer,
            majority: majority_cutoff,
            end_join,
            end_join_option,
            overlap: overlap_size,
            tcs_qc,
            ref_genome,
            ref_start,
            ref_start_lower,
            ref_end,
            ref_end_lower,
            indel,
            trim,
            trim_ref,
            trim_ref_start,
            trim_ref_end,
        });

        print!("Add another region? (y/n, default as n):\n>  ");
        match collect_input().as_str() {
            "y" | "Y" => continue,
            _ => break,
        };
    }

    print!("Enter your email address (optional):\n>  ");
    let email = collect_input();
    let email = if email.is_empty() { None } else { Some(email) };

    let params = Params {
        platform_error_rate: error_rate,
        platform_format: platform,
        email: email,
        primer_pairs: regions,
    };

    println!("Your input directory: {}", input_dir);
    println!("Your entered parameters: ");
    println!("{}", params);

    print!("\nDo you wish to save the parameters to a JSON file? (y/n):\n>  ");
    let save = match collect_input().as_str() {
        "y" | "Y" => true,
        _ => false,
    };

    if save {
        let json = serde_json::to_string_pretty(&params).expect("Failed to serialize");
        loop {
            print!("Enter the path to save the JSON file (e.g. /path/to/params.json):\n>  ");
            let json_path = PathBuf::from(collect_input());

            let mut file = match File::create(&json_path) {
                Ok(file) => {
                    println!("File created successfully at {}.", json_path.display());
                    file
                }
                Err(e) => {
                    eprintln!("Error creating file: {}, retry", e);
                    continue;
                }
            };

            match file.write_all(json.as_bytes()) {
                Ok(_) => {
                    println!("Parameters saved to JSON file at {}.", json_path.display());
                    break;
                }
                Err(e) => {
                    eprintln!("Error writing to file: {}, retry", e);
                    continue;
                }
            }
        }
    } else {
        println!("Parameters not saved. Goodbye!");
    }
}

fn collect_input() -> String {
    io::stdout().flush().unwrap();
    let mut input = String::new();
    io::stdin()
        .read_line(&mut input)
        .expect("Failed to read line");
    input.trim().to_string()
}

fn collect_end_join_option() -> u8 {
    loop {
        print!(
            "End-join option? Choose from (1-4):\n\
            1: simple join, no overlap\n\
            2: known overlap\n\
            3: unknow overlap, use sample consensus to determine overlap, all sequence pairs have same overlap\n\
            4: unknow overlap, determine overlap by individual sequence pairs, sequence pairs can have different overlap\n\
            >  "
        );
        match collect_input().parse::<u8>() {
            Ok(num) if (1..=4).contains(&num) => break num,
            _ => {
                eprintln!("Invalid input. Please enter a number from 1 to 4.");
                continue;
            }
        }
    }
}

fn get_ref() -> String {
    loop {
        print!(
            "Choose reference genome (1-3):\n\
        1: HIV-1 HXB2\n\
        2: SIV MAC239\n\
        >  "
        );
        match collect_input().parse::<u8>() {
            Ok(num) if (1..=2).contains(&num) => {
                if num == 1 {
                    break "HXB2".to_string();
                } else {
                    break "SIVmm239".to_string();
                }
            }
            _ => {
                eprintln!("Invalid input. Please enter a number from 1 to 2.");
                continue;
            }
        }
    }
}

fn get_ref_and_locations() -> (String, u32, Option<u32>, u32, Option<u32>) {
    let ref_genome = get_ref();

    print!("reference 5'end ref position range starts, 0 if no need to match this end \n>  ");
    let ref_start = match collect_input().as_str() {
        "" => 0,
        input => input.parse::<u32>().unwrap_or(0),
    };

    print!("reference 5'end ref position range ends, 0 if no need to match this end \n>  ");
    let ref_start_lower = match collect_input().as_str() {
        "" => 0,
        input => input.parse::<u32>().unwrap_or(0),
    };

    print!("reference 3'end ref position range starts: 0 if no need to match this end \n>  ");
    let ref_end = match collect_input().as_str() {
        "" => 0,
        input => input.parse::<u32>().unwrap_or(0),
    };

    print!("reference 3'end ref position range ends: 0 if no need to match this end \n>  ");
    let ref_end_lower = match collect_input().as_str() {
        "" => 0,
        input => input.parse::<u32>().unwrap_or(0),
    };
    (
        ref_genome,
        ref_start,
        Some(ref_start_lower),
        ref_end,
        Some(ref_end_lower),
    )
}
