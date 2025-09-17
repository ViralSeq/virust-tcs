#![allow(non_snake_case)]

use std::collections::HashMap;
use std::error::Error;

use getset::{Getters, Setters};
use serde::{Deserialize, Serialize};

use crate::helper::drm_helper::{Coord, DrmList, DrmRefInfo, DrmVersion};

// This is the structure for DRM region configuration used in the SDRM pipeline
// It is populated based on the selected DRM version (DrmVersion), and the DrmList
// This struct is then used in the SDRM pipeline to interpret sequence data and identify DRMs

#[derive(Debug, PartialEq, Eq, Getters, Setters, Serialize, Deserialize)]
pub struct DrmRegionConfig {
    #[getset(get = "pub", set = "pub")]
    drm_version: String,
    #[getset(get = "pub", set = "pub")]
    region: String,
    #[getset(get = "pub", set = "pub")]
    drm_classes: Vec<String>,
    #[getset(get = "pub", set = "pub")]
    drm_classes_with_range: HashMap<String, Vec<u32>>,
    #[getset(get = "pub", set = "pub")]
    drm_list: DrmList,
    #[getset(get = "pub", set = "pub")]
    seq_coord: Coord,
    #[getset(get = "pub", set = "pub")]
    ref_info: DrmRefInfo,
}

impl DrmRegionConfig {
    pub fn from_drm_version(
        drm_version: &DrmVersion,
        drm_master_list: &DrmList,
        region_name: &str,
    ) -> Result<Self, Box<dyn Error>> {
        let drm_classes: Vec<String> = drm_version
            .seq_drm_correlation()
            .get(region_name)
            .cloned()
            .ok_or(format!(
                "Region name {} not found in seq_drm_correlation",
                region_name
            ))?;

        let mut drm_classes_with_range: HashMap<String, Vec<u32>> = HashMap::new();

        let mut drm_list = DrmList::new();

        for drm_class in &drm_classes {
            let drm_range = drm_version
                .DRM_range()
                .get(drm_class.as_str())
                .ok_or(format!("DRM class {} not found in DRM_range", drm_class))?;

            drm_classes_with_range.insert(drm_class.clone(), drm_range.clone());

            let mut drm_single_class_list = drm_master_list.get(drm_class).cloned().ok_or(
                format!("DRM class {} not found in master DRM list", drm_class),
            )?;

            drm_single_class_list.retain(|m| drm_range.contains(&m.position()));
            drm_list.insert(drm_class.clone(), drm_single_class_list);
        }

        let seq_coord = drm_version
            .seq_coord()
            .get(region_name)
            .cloned()
            .ok_or(format!(
                "Region name {} not found in seq_coord",
                region_name
            ))?;

        let ref_info = drm_version
            .ref_info()
            .get_one_region(region_name)
            .ok_or(format!("Region name {} not found in ref_info", region_name))?;

        let drm_region_config = DrmRegionConfig {
            drm_version: drm_version.version().clone(),
            region: region_name.to_string(),
            drm_classes,
            drm_classes_with_range,
            drm_list: drm_list,
            seq_coord: seq_coord,
            ref_info: ref_info,
        };

        Ok(drm_region_config)
    }
}

#[cfg(test)]
mod tests {

    use super::*;
    use crate::helper::drm_helper::DrmListTrait;

    #[test]
    fn test_drm_region_config_from_drm_version() {
        let drm_version = DrmVersion::build_from_version("v4").unwrap();
        let drm_master_list = DrmList::build().unwrap();
        let region_name = "IN";
        let drm_region_config =
            DrmRegionConfig::from_drm_version(&drm_version, &drm_master_list, region_name);
        assert!(drm_region_config.is_ok());
        let drm_region_config = drm_region_config.unwrap();
        dbg!(&drm_region_config);
    }
}
