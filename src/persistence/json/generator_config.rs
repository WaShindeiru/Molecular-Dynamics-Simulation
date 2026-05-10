use std::{fs, io};

use crate::data::ValueUnits;
use crate::simulations::generators::generator_config::dense::DenseGeneratorConfig;
use crate::simulations::generators::generator_config::GeneratorConfig;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct GeneratorConfigFile {
  #[serde(default)]
  pub value_units: ValueUnits,
  pub generator: GeneratorConfig,
}

impl GeneratorConfigFile {
  pub fn new(generator: GeneratorConfig, value_units: ValueUnits) -> Self {
    Self {
      value_units,
      generator,
    }
  }

  pub fn to_value_units(&self, target: ValueUnits) -> Self {
    Self {
      value_units: target,
      generator: self.generator.to_value_units(self.value_units, target),
    }
  }

  pub fn into_generator_config_unitless(self) -> GeneratorConfig {
    self
      .to_value_units(ValueUnits::Unitless)
      .generator
  }

  pub fn to_json_string(&self) -> Result<String, serde_json::Error> {
    serde_json::to_string_pretty(self)
  }

  pub fn from_json_str(s: &str) -> Result<Self, serde_json::Error> {
    match serde_json::from_str::<Self>(s) {
      Ok(value) => Ok(value),
      Err(_) => {
        // Backward compatibility: old format with only GeneratorConfig.
        let generator = serde_json::from_str::<GeneratorConfig>(s)?;
        Ok(Self {
          value_units: ValueUnits::Unitless,
          generator,
        })
      }
    }
  }

  pub fn to_json_file(&self, path: &str) -> io::Result<()> {
    let json = self
      .to_json_string()
      .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;
    fs::write(path, json)
  }

  pub fn from_json_file(path: &str) -> io::Result<Self> {
    let content = fs::read_to_string(path)?;
    Self::from_json_str(&content).map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))
  }
}
