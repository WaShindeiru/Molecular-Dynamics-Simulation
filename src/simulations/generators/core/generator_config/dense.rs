use nalgebra::Vector3;
use std::{fs, io};

use crate::data::SimulationConfig;
use crate::data::units::{R_U, ValueUnits};
use crate::simulations::generators::core::generate::dense::DenseGenerator;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct DenseGeneratorConfig {
  pub particle_distance: f64,
  pub offset: Vector3<f64>,
}

impl DenseGeneratorConfig {
  pub fn new(particle_distance: f64, offset: Vector3<f64>) -> Self {
    Self {
      particle_distance,
      offset,
    }
  }

  pub fn to_generator(&self, simulation_config: &SimulationConfig) -> DenseGenerator {
    DenseGenerator::new(
      simulation_config.potential_gravity_max,
      simulation_config.world_size,
      self.particle_distance,
      self.offset,
    )
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    let scale = ValueUnits::scale_between(source, target, R_U);

    Self {
      particle_distance: self.particle_distance * scale,
      offset: self.offset * scale,
    }
  }
}
