use crate::data::SimulationConfig;
use crate::data::units::{VELOCITY_U, ValueUnits};
use crate::persistence::json::particle_config::Vector3Record;
use crate::simulations::generators::core::generate::random::RandomGenerator;

#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RandomGeneratorConfig {
  pub num_particles: usize,
  pub vel_mean: Vector3Record,
  pub vel_std_dev: Vector3Record,
}

impl RandomGeneratorConfig {
  pub fn new(
    num_particles: usize,
    vel_mean: Vector3Record,
    vel_std_dev: Vector3Record,
  ) -> Self {
    Self {
      num_particles,
      vel_mean,
      vel_std_dev,
    }
  }

  pub fn to_generator(&self, simulation_config: &SimulationConfig) -> RandomGenerator {
    RandomGenerator::new(
      self.num_particles,
      simulation_config.initial_gravity(),
      simulation_config.world_size,
      self.vel_mean.to_runtime(),
      self.vel_std_dev.to_runtime(),
    )
  }

  pub fn to_value_units(&self, source: ValueUnits, target: ValueUnits) -> Self {
    let vel_scale = ValueUnits::scale_between(source, target, VELOCITY_U);

    Self {
      num_particles: self.num_particles,
      vel_mean: self.vel_mean.scale(vel_scale),
      vel_std_dev: self.vel_std_dev.scale(vel_scale),
    }
  }
}
