#![allow(non_snake_case)]

use std::collections::HashMap;
use std::error::Error;

use getset::{Getters, Setters};
use itertools::Itertools;
use serde::{Deserialize, Serialize};

// DRM version configuration structures and functions
// DRM versions can be edited in the JSON file at resources/drm_config/drm_versions_config.json
// The structure here should match the JSON structure
// Add new DRM versions by editing the JSON file, no need to change the code
// DrmVersion is corresponding to each library prepration kit version, which covers different regions of HIV genome, and potentially different DRM positions.
#[derive(Debug, PartialEq, Eq, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct DrmVersion {
    #[getset(get = "pub", set = "pub")]
    version: String,
    #[getset(get = "pub", set = "pub")]
    DRM_range: DRMRange,
    #[getset(get = "pub", set = "pub")]
    seq_coord: SeqCoord,
    #[getset(get = "pub", set = "pub")]
    seq_drm_correlation: HashMap<String, Vec<String>>,
    #[getset(get = "pub", set = "pub")]
    ref_info: DrmRefInfo,
}

impl DrmVersion {
    pub fn build_from_version(version: &str) -> Result<Self, Box<dyn Error>> {
        get_drm_version(version)
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct DRMRange {
    #[getset(get = "pub", set = "pub")]
    CAI: Vec<u32>,
    #[getset(get = "pub", set = "pub")]
    PI: Vec<u32>,
    #[getset(get = "pub", set = "pub")]
    NRTI: Vec<u32>,
    #[getset(get = "pub", set = "pub")]
    NNRTI: Vec<u32>,
    #[getset(get = "pub", set = "pub")]
    INSTI: Vec<u32>,
}

impl DRMRange {
    pub fn get(&self, key: &str) -> Option<&Vec<u32>> {
        match key {
            "CAI" => Some(&self.CAI),
            "PI" => Some(&self.PI),
            "NRTI" => Some(&self.NRTI),
            "NNRTI" => Some(&self.NNRTI),
            "INSTI" => Some(&self.INSTI),
            _ => None,
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct SeqCoord {
    #[getset(get = "pub", set = "pub")]
    CA: Coord,
    #[getset(get = "pub", set = "pub")]
    PR: Coord,
    #[getset(get = "pub", set = "pub")]
    RT: Coord,
    #[getset(get = "pub", set = "pub")]
    IN: Coord,
}

impl SeqCoord {
    pub fn get(&self, key: &str) -> Option<&Coord> {
        match key {
            "CA" => Some(&self.CA),
            "PR" => Some(&self.PR),
            "RT" => Some(&self.RT),
            "IN" => Some(&self.IN),
            _ => None,
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct Coord {
    #[getset(get = "pub", set = "pub")]
    minimum: u32,
    #[getset(get = "pub", set = "pub")]
    maximum: u32,
    #[getset(get = "pub", set = "pub")]
    gap: Option<Gap>,
}

#[derive(Debug, PartialEq, Eq, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct DrmRefInfo {
    #[getset(get = "pub", set = "pub")]
    ref_type: String,
    #[getset(get = "pub", set = "pub")]
    ref_coord: HashMap<String, [u32; 2]>,
}

impl DrmRefInfo {
    pub fn get_one_region(&self, region_name: &str) -> Option<Self> {
        if let Some(coord) = self.ref_coord.get(region_name) {
            Some(DrmRefInfo {
                ref_type: self.ref_type.clone(),
                ref_coord: HashMap::from([(region_name.to_string(), *coord)]),
            })
        } else {
            None
        }
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct Gap {
    #[getset(get = "pub", set = "pub")]
    minimum: u32,
    #[getset(get = "pub", set = "pub")]
    maximum: u32,
}

pub fn get_drm_version_config() -> Result<HashMap<String, DrmVersion>, Box<dyn Error>> {
    let drm_version_json =
        include_str!("../../../resources/drm_config/drm_versions_config.json").to_string();
    let configs: Vec<DrmVersion> = serde_json::from_str(&drm_version_json)?;

    let config_map: HashMap<String, DrmVersion> = configs
        .into_iter()
        .map(|cfg| (cfg.version.clone(), cfg))
        .collect();
    Ok(config_map)
}

pub fn get_supported_drm_versions() -> Result<Vec<String>, Box<dyn Error>> {
    let config_map = get_drm_version_config()?;
    let versions: Vec<String> = config_map.keys().cloned().collect();
    Ok(versions)
}

pub fn get_drm_version(version: &str) -> Result<DrmVersion, Box<dyn Error>> {
    let config_map = get_drm_version_config()?;
    let mut version = version.to_lowercase();
    if version == "v2" {
        version = "v1".to_string(); // this is because v2 is identical to v1 in terms of DRM config
    }

    match config_map.get(&version) {
        Some(v) => Ok(v.clone()),
        None => Err(format!(
            "DRM version not found. Supported versions are: { }. Note that 'v2' is an alias for 'v1'.",
            config_map
                .keys()
                .cloned()
                .sorted()
                .collect::<Vec<String>>()
                .join(", ")
        )
        .into()),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_drm_version_config() {
        let config = get_drm_version_config();
        dbg!(config.as_ref().unwrap());
        assert!(config.is_ok());
    }

    #[test]
    fn test_get_supported_drm_versions() {
        let versions = get_supported_drm_versions();
        dbg!(versions.as_ref().unwrap().iter().sorted());
        assert!(versions.is_ok());
        assert!(versions.unwrap().contains(&"v1".to_string()));
    }

    #[test]
    fn test_get_drm_version() {
        let drm_v1 = get_drm_version("v1");
        assert!(drm_v1.is_ok());
        let drm_v1 = drm_v1.unwrap();
        dbg!(&drm_v1);
        assert_eq!(drm_v1.version, "v1".to_string());

        dbg!(&drm_v1.ref_info().get_one_region("CA"));
        assert!(drm_v1.ref_info().get_one_region("CA").is_some());
        assert!(drm_v1.ref_info().get_one_region("XX").is_none());

        let drm_non_exist = get_drm_version("vx");
        dbg!(drm_non_exist.as_ref().err());
        assert!(drm_non_exist.is_err());
    }
}
