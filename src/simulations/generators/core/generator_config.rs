pub mod dense;
pub mod nanotube;
pub mod velocity_nanotube;

use crate::data::SimulationConfig;
use crate::data::ValueUnits;
use crate::simulations::generators::core::generate::GeneratorType;
use crate::simulations::generators::core::generator_config::dense::DenseGeneratorConfig;
use crate::simulations::generators::core::generator_config::nanotube::NanotubeGeneratorConfig;
use crate::simulations::generators::core::generator_config::velocity_nanotube::VelocityNanotubeGeneratorConfig;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
#[serde(tag = "type", content = "config", rename_all = "snake_case")]
pub enum GeneratorConfig {
  Dense(DenseGeneratorConfig),
  Nanotube(NanotubeGeneratorConfig),
  VelocityNanotube(VelocityNanotubeGeneratorConfig),
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
      GeneratorConfig::VelocityNanotube(config) => {
        GeneratorType::VelocityNanotube(config.to_generator())
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
      GeneratorConfig::VelocityNanotube(config) => {
        GeneratorConfig::VelocityNanotube(config.to_value_units(source, target))
      }
    }
  }
}
