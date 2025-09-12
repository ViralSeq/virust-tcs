//TODO: Implement FastQC analysis functions

use std::error::Error;
use std::fs::File;
use std::io::BufReader;
use std::path::Path;

use bio::io::fastq;
use getset::{self, Getters, Setters};
use plotters::prelude::*;
use serde::{Deserialize, Serialize};
use statrs::statistics::{Data, Distribution, Max, Min, OrderStatistics};

#[derive(Debug, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct FastQcResults {
    #[getset(get = "pub", set = "pub")]
    total_reads: usize,
    #[getset(get = "pub", set = "pub")]
    read_length: Vec<usize>,
    #[getset(get = "pub", set = "pub")]
    quality_score_distribution: Vec<QualityScoreDistribution>,
}

#[derive(Debug, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct QualityScoreDistribution {
    // notice this position starts from 1, not 0.
    #[getset(get = "pub", set = "pub")]
    position: usize,
    #[getset(get = "pub", set = "pub")]
    count: usize,
    #[getset(get = "pub", set = "pub")]
    quality_mean: f64,
    #[getset(get = "pub", set = "pub")]
    quality_min: f64,
    #[getset(get = "pub", set = "pub")]
    quality_max: f64,
    #[getset(get = "pub", set = "pub")]
    quality_first_quartile: f64,
    #[getset(get = "pub", set = "pub")]
    quality_median: f64,
    #[getset(get = "pub", set = "pub")]
    quality_third_quartile: f64,
    #[getset(get = "pub", set = "pub")]
    quality_standard_deviation: f64,
}

impl FastQcResults {
    pub fn new() -> Self {
        FastQcResults {
            total_reads: 0,
            read_length: Vec::new(),
            quality_score_distribution: Vec::new(),
        }
    }

    pub fn export_quality_score_distribution_to_csv(
        &self,
        path_to_csv: &Path,
    ) -> Result<(), Box<dyn Error>> {
        let mut wtr = csv::Writer::from_path(path_to_csv)?;

        for (_i, qsd) in self.quality_score_distribution.iter().enumerate() {
            wtr.serialize(qsd)?;
        }

        wtr.flush()?;
        Ok(())
    }
}

impl QualityScoreDistribution {
    pub fn new() -> Self {
        QualityScoreDistribution {
            position: 0,
            count: 0,
            quality_mean: 0.0,
            quality_min: 0.0,
            quality_max: 0.0,
            quality_first_quartile: 0.0,
            quality_median: 0.0,
            quality_third_quartile: 0.0,
            quality_standard_deviation: 0.0,
        }
    }

    pub fn from_qual_vec(qual_vec: Vec<u8>) -> Self {
        let mut data = Data::new(qual_vec.iter().map(|&q| q as f64).collect::<Vec<f64>>());

        let quality_mean = data.mean().unwrap_or(0.0);

        let quality_min = data.min();
        let quality_max = data.max();
        // notice that the quantile estimation is using Type 8 in the `quantile` function in R.
        // The resulting quantile estimates are approximately median-unbiased regardless of the distribution of x.
        let quality_first_quartile = data.quantile(0.25);
        let quality_median = data.median();
        let quality_third_quartile = data.quantile(0.75);
        let quality_standard_deviation = data.std_dev().unwrap_or(0.0);

        QualityScoreDistribution {
            position: 0,
            count: qual_vec.len(),
            quality_mean,
            quality_min,
            quality_max,
            quality_first_quartile,
            quality_median,
            quality_third_quartile,
            quality_standard_deviation,
        }
    }
}

pub fn fastqc_analysis(fastq_file_path: &Path) -> Result<FastQcResults, Box<dyn Error>> {
    let mut results = FastQcResults::new();

    let file = File::open(fastq_file_path)?;
    let reader = fastq::Reader::new(BufReader::new(file));

    let mut qual_scores: Vec<Vec<u8>> = Vec::new();

    for record in reader.records() {
        let record = record?;
        let qual = record
            .qual()
            .iter()
            .map(|q| q.saturating_sub(33))
            .collect::<Vec<u8>>(); // Important: Convert ASCII to Phred quality scores
        qual_scores.push(qual);
    }

    results.total_reads = qual_scores.len();

    let length_distribution = qual_scores.iter().map(|q| q.len()).collect::<Vec<usize>>();
    results.read_length = length_distribution.clone();

    let max_length = length_distribution.iter().max().cloned().unwrap_or(0);

    for i in 0..max_length {
        let mut qual_vec = Vec::new();
        for q in qual_scores.iter() {
            if i < q.len() {
                let adjusted_qual_score = q[i].clamp(0, 40); // Max at 40, Q40 allowed

                qual_vec.push(adjusted_qual_score);
            }
        }

        let mut qds = QualityScoreDistribution::from_qual_vec(qual_vec);

        qds.position = i + 1;
        results.quality_score_distribution.push(qds);
    }

    Ok(results)
}

// function to plot the quality score distribution using plotters crate

pub fn plot_quality_score_distribution(
    distribution: &[QualityScoreDistribution],
    output_path: &Path,
) -> Result<(), Box<dyn Error>> {
    let root = BitMapBackend::new(output_path, (1200, 800)).into_drawing_area();
    root.fill(&WHITE)?;

    let mut chart = ChartBuilder::on(&root)
        .caption("Quality Score Distribution", ("sans-serif", 50).into_font())
        .margin(10)
        .x_label_area_size(50)
        .y_label_area_size(50)
        .build_cartesian_2d(1f64..(distribution.len() as f64), 0f64..40f64)?;

    // Configure mesh (grid lines) and axis labels
    chart
        .configure_mesh()
        .x_desc("Position in Read")
        .y_desc("Quality Score")
        .axis_desc_style(("sans-serif", 20))
        .x_label_formatter(&|x| format!("{}", x))
        .y_label_formatter(&|y| format!("{:.0}", y))
        .draw()?;

    // Draw line for median values
    chart.draw_series(LineSeries::new(
        distribution
            .iter()
            .enumerate()
            .map(|(i, qsd)| ((i as f64), qsd.quality_median)),
        BLUE.stroke_width(2),
    ))?;

    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_quality_score_distribution_from_vec() {
        let qual_vec = vec![30, 32, 28, 35, 31];
        let qsd = QualityScoreDistribution::from_qual_vec(qual_vec);
        assert_eq!(qsd.quality_mean, 31.2);
        assert_eq!(qsd.quality_min, 28.0);
        assert_eq!(qsd.quality_max, 35.0);
        assert_eq!((qsd.quality_first_quartile * 10.0).round() / 10.0, 29.3);
        assert_eq!(qsd.quality_median, 31.0);
        assert_eq!(qsd.quality_third_quartile, 33.0);
        assert!((qsd.quality_standard_deviation - 2.588436).abs() < 0.001); // Approximate value
    }

    #[test]
    fn test_fastqc_analysis() {
        let fastq_file_path = Path::new("tests/data/test_fastqc/sample.fastq");
        let results = fastqc_analysis(&fastq_file_path);
        assert!(results.is_ok());
        let results = results.unwrap();
        assert_eq!(results.total_reads, 907);
        assert_eq!(results.read_length.len(), 907);
        assert_eq!(results.quality_score_distribution.len(), 524);

        //write tests for the quality score distribution to json file
        let json = serde_json::to_string_pretty(&results).unwrap();
        std::fs::write("tests/data/test_fastqc/sample.json", json).unwrap();

        results
            .export_quality_score_distribution_to_csv(Path::new(
                "tests/data/test_fastqc/sample.csv",
            ))
            .unwrap();

        plot_quality_score_distribution(
            &results.quality_score_distribution(),
            Path::new("tests/data/test_fastqc/sample.png"),
        )
        .unwrap()
    }
}
