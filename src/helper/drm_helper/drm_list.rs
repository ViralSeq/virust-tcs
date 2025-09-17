#![allow(non_snake_case)]

use std::collections::HashMap;
use std::error::Error;

use getset::{Getters, Setters};
use serde::{Deserialize, Serialize};

// DRM list structure and function
// The DRM list can be edited in the JSON file at resources/drm_config/drm_list.json
// The structure here should match the JSON structure
// Add new DRM mutations by editing the JSON file, no need to change the code
// HIV DRM list is adapted from Stanford HIVdb (https://hivdb.stanford.edu/)
// Future adaption to other versions of HIVdb or other viruses can be done by editing the JSON file
// DRMs to CAIs may neeed to be updated frequently as new drugs are approved and new DRMs are identified
// Future versions of the DRM List should be able to use HIVdb versioning system if available

pub type DrmList = HashMap<String, Vec<Mutation>>;

pub trait DrmListTrait {
    fn get_classes(&self, class: &str) -> Option<&[Mutation]>;
    fn find(&self, class: &str, position: u32) -> Option<&Mutation>;
    fn build() -> Result<DrmList, Box<dyn Error>>;
}

impl DrmListTrait for DrmList {
    fn get_classes(&self, class: &str) -> Option<&[Mutation]> {
        self.get(class).map(|v| v.as_slice())
    }

    fn find(&self, class: &str, position: u32) -> Option<&Mutation> {
        self.get(class)
            .and_then(|mutations| mutations.iter().find(|m| m.position == position))
    }

    fn build() -> Result<DrmList, Box<dyn Error>> {
        let drm_list_str = include_str!("../../../resources/drm_config/drm_list.json");
        let drm_list: DrmList = serde_json::from_str(drm_list_str)?;
        Ok(drm_list)
    }
}

#[derive(Debug, PartialEq, Eq, Clone, Getters, Setters, Serialize, Deserialize)]
pub struct Mutation {
    #[getset(get = "pub", set = "pub")]
    position: u32,
    #[getset(get = "pub", set = "pub")]
    #[serde(rename = "wild-type")]
    wild_type: String,
    #[getset(get = "pub", set = "pub")]
    mutations: Vec<String>,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_get_drm_list() {
        let drm_list = DrmList::build();

        assert!(drm_list.is_ok());
        assert!(drm_list.unwrap().contains_key("CAI"));
    }
}
