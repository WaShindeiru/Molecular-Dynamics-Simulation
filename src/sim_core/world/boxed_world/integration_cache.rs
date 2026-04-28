pub mod integration_cache_builder;

use std::collections::HashMap;
use std::sync::Arc;
use nalgebra::Vector3;
use crate::sim_core::world::boundary_constraint::ParticleCompliance;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;

pub struct IntegrationCache {
  pub integration_box_cache: BoxContainer<Arc<SimulationBox>>,
  pub integration_half_velocity_cache: HashMap<usize, Vector3<f64>>,
  pub particle_compliance: HashMap<usize, ParticleCompliance>,
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
}