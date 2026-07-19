pub mod integration_cache_builder;

use std::sync::Arc;

use nalgebra::Vector3;

use crate::sim_core::world::boundary_constraint::ParticleCompliance;
use crate::sim_core::world::cell::LinkedCellContainer;

pub struct IntegrationCache {
  pub local_container: LinkedCellContainer,
  pub read_container: Arc<LinkedCellContainer>,
  pub half_velocity_cache: Vec<Vector3<f64>>,
  pub particle_compliance: Vec<ParticleCompliance>,
}

impl IntegrationCache {
  pub fn new(
    local_container: LinkedCellContainer,
    read_container: Arc<LinkedCellContainer>,
    half_velocity_cache: Vec<Vector3<f64>>,
    particle_compliance: Vec<ParticleCompliance>,
  ) -> Self {
    IntegrationCache { local_container, read_container, half_velocity_cache, particle_compliance }
  }
}
