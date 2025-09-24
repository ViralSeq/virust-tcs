use serde::{ Deserialize, Serialize };
use std::collections::BTreeMap;

/// Top-level: many libraries keyed by ID (e.g., "RV95")
pub type Libraries = BTreeMap<String, LibData>;

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LibData {
    pub raw_distribution: RawDistribution, // Vec<(String, u32)>
    pub raw_sequence_analysis: RawSequenceAnalysis, // Vec<(String, u32)> + drilldowns
    pub number_at_regions: NumberAtRegions, // Vec<(String, u32, u32, u32)>
    pub distinct_to_raw: RegionF32Table, // Vec<(String, f32)>
    pub resampling_index: RegionF32Table, // Vec<(String, f32)>
    pub size_distribution: SizeDistribution, // Map<String, Vec<(u32, Option<f32>, Option<u32>)>>
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawDistribution {
    pub data: Vec<(String, u32)>, // (region, paired_raw)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RawSequenceAnalysis {
    pub data: Vec<(String, u32)>, // (category, count)
    pub drilldowns: Vec<DrilldownU32>, // label + same (String, u32) rows
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct DrilldownU32 {
    pub label: String,
    pub data: Vec<(String, u32)>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NumberAtRegions {
    pub data: Vec<(String, Option<u32>, Option<u32>, Option<u32>)>, // (region, tcs, combined_tcs, tcs_after_qc)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RegionF32Table {
    pub data: Vec<(String, f32)>, // (region, value)
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SizeDistribution {
    pub data: BTreeMap<String, Vec<(u32, Option<f32>, Option<u32>)>>, // region -> Vec<(index, distribution, cutoff)>
}
