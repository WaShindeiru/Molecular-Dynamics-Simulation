pub mod integration_cache_builder;

use crate::sim_core::world::boundary_constraint::ParticleCompliance;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;

pub struct IntegrationCache {
  integration_box_cache: BoxContainer<Arc<SimulationBox>>,
  integration_half_velocity_cache: HashMap<usize, Vector3<f64>>,
  particle_compliance: HashMap<usize, ParticleCompliance>,
}

impl IntegrationCache {
  pub fn new(
    integration_box_cache: BoxContainer<Arc<SimulationBox>>,
    integration_half_velocity_cache: HashMap<usize, Vector3<f64>>,
    particle_compliance: HashMap<usize, ParticleCompliance>,
  ) -> Self {
    IntegrationCache {
      integration_box_cache,
      integration_half_velocity_cache,
      particle_compliance,
    }
  }

  pub fn box_cache(&self) -> &BoxContainer<Arc<SimulationBox>> {
    &self.integration_box_cache
  }

  pub fn half_velocity_cache(&self) -> &HashMap<usize, Vector3<f64>> {
    &self.integration_half_velocity_cache
  }

  pub fn particle_compliance(&self) -> &HashMap<usize, ParticleCompliance> {
    &self.particle_compliance
  }
}
