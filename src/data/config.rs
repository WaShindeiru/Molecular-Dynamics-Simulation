pub mod builder;
pub mod particle_config;
pub mod simulation_config;

use crate::data::SimulationConfig;
use particle_config::ParticleConfig;

#[derive(Clone)]
pub struct ConfigAll {
  pub simulation_config: SimulationConfig,
  pub particle_config: ParticleConfig,
}

impl ConfigAll {
  pub fn new(simulation_config: SimulationConfig, particle_config: ParticleConfig) -> Self {
    Self {
      simulation_config,
      particle_config,
    }
  }
}
