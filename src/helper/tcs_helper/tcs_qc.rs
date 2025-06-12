use std::collections::HashMap;
use std::error::Error;

use getset::{Getters, Setters};
use virust_locator::config::Args as LocatorArgs;
use virust_locator::locator::Locator;

#[derive(Debug, Clone, Getters, Setters)]
pub struct TcsQcInput {
    #[getset(get = "pub", set = "pub")]
    query: Vec<String>,
    #[getset(get = "pub", set = "pub")]
    reference: QcReference,
    #[getset(get = "pub", set = "pub")]
    algorithm: QcAlgorithm,
}

impl TcsQcInput {
    pub fn with_attrs(query: Vec<&[u8]>, reference: String, algorithm_code: u8) -> Option<Self> {
        let reference = QcReference::from_string(&reference).unwrap_or_default();
        let algorithm = QcAlgorithm::from_option_code(algorithm_code).unwrap_or_default();
        if query.is_empty() {
            return None;
        }
        let query = query
            .into_iter()
            .map(|q| String::from_utf8_lossy(q).to_string())
            .collect();
        Some(TcsQcInput {
            query,
            reference,
            algorithm,
        })
    }

    pub fn to_locator_args(&self) -> LocatorArgs {
        LocatorArgs {
            query: self.query.clone(),
            reference: self.reference.to_string(),
            type_query: "nt".to_string(),
            algorithm: self.algorithm.to_option_code(),
        }
    }

    //TODO: Implement the actual locator run logic
    pub fn run_locator(&self) -> Result<TcsQcOutput, Box<dyn Error + Send + Sync>> {
        let args = self.to_locator_args();
        let locator = Locator::build(&args)?;
        let mut query_locator_hashmap: HashMap<&[u8], Option<Locator>> = HashMap::new();
        for i in 0..self.query().len() {
            let query = self.query()[i].as_bytes();
            query_locator_hashmap.insert(query, locator[i].clone());
        }
        let tcs_qc_output = TcsQcOutput {
            results_map: query_locator_hashmap,
        };
        Ok(tcs_qc_output)
    }
}

#[derive(Debug, Clone)]
pub enum QcAlgorithm {
    SemiGlobal,
    PatternMatching,
}
impl Default for QcAlgorithm {
    fn default() -> Self {
        QcAlgorithm::SemiGlobal
    }
}

impl QcAlgorithm {
    pub fn to_option_code(&self) -> u8 {
        match self {
            QcAlgorithm::SemiGlobal => 1,
            QcAlgorithm::PatternMatching => 2,
        }
    }

    pub fn from_option_code(code: u8) -> Option<Self> {
        match code {
            1 => Some(QcAlgorithm::SemiGlobal),
            2 => Some(QcAlgorithm::PatternMatching),
            _ => None,
        }
    }
}

#[derive(Debug, Clone)]
pub enum QcReference {
    HXB2,
    SIVmm239,
}

impl Default for QcReference {
    fn default() -> Self {
        QcReference::HXB2
    }
}

impl QcReference {
    pub fn to_string(&self) -> String {
        match self {
            QcReference::HXB2 => "HXB2".to_string(),
            QcReference::SIVmm239 => "SIVmm239".to_string(),
        }
    }

    pub fn from_string(reference: &str) -> Option<Self> {
        match reference {
            "HXB2" => Some(QcReference::HXB2),
            "SIVmm239" => Some(QcReference::SIVmm239),
            _ => None,
        }
    }
}

#[derive(Debug, Clone, Getters, Setters)]
pub struct TcsQcOutput<'a> {
    #[getset(get = "pub")]
    results_map: HashMap<&'a [u8], Option<Locator>>,
}
