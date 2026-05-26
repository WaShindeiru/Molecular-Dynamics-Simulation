pub mod dense;
pub mod nanotube;

use crate::data::SimulationConfig;
use crate::data::ValueUnits;
use crate::simulations::generators::generate::GeneratorType;
use crate::simulations::generators::generator_config::dense::DenseGeneratorConfig;
use crate::simulations::generators::generator_config::nanotube::NanotubeGeneratorConfig;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type", content = "config", rename_all = "snake_case")]
pub enum GeneratorConfig {
  Dense(DenseGeneratorConfig),
  Nanotube(NanotubeGeneratorConfig),
}

impl GeneratorConfig {
  pub fn to_generator(&self, simulation_config: &SimulationConfig) -> GeneratorType {
    match self {
      GeneratorConfig::Dense(config) => {
        GeneratorType::Dense(config.to_generator(simulation_config))
      }
      GeneratorConfig::Nanotube(config) => {
        GeneratorType::Nanotube(config.to_generator(simulation_config))
      }
    }
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    match self {
      GeneratorConfig::Dense(config) => {
        GeneratorConfig::Dense(config.to_value_units(source, target))
      }
      GeneratorConfig::Nanotube(config) => {
        GeneratorConfig::Nanotube(config.to_value_units(source, target))
      }
    }
  }
}
