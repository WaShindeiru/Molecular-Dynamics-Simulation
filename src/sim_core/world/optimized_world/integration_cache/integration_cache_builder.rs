use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::sim_core::world::boundary_constraint::{Compliance, ParticleCompliance};
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::cell::LinkedCellContainer;
use crate::sim_core::world::optimized_world::integration_cache::IntegrationCache;

/// Never read before being overwritten: `add_velocity_results` always writes the real value
/// for id `i` in the same loop iteration where the backing container's `present[i]` becomes
/// `true`, and `build()` only exposes these vectors once both containers report `is_complete`.
fn placeholder_compliance() -> ParticleCompliance {
  ParticleCompliance {
    compliant: false,
    x: Compliance::Compliant,
    y: Compliance::Compliant,
    z: Compliance::Compliant,
  }
}

pub struct IntegrationCacheBuilder {
  old: Arc<LinkedCellContainer>,
  local_container: LinkedCellContainer,
  read_container: LinkedCellContainer,
  half_velocity: Vec<Vector3<f64>>,
  particle_compliance: Vec<ParticleCompliance>,
}

impl IntegrationCacheBuilder {
  pub fn new(old: Arc<LinkedCellContainer>) -> Self {
    let num_particles = old.particles().len();
    let config = *old.config();
    let edge_condition = old.edge_condition();
    IntegrationCacheBuilder {
      local_container: LinkedCellContainer::new_empty(num_particles, config, edge_condition),
      read_container: LinkedCellContainer::new_empty(num_particles, config, edge_condition),
      half_velocity: vec![Vector3::zeros(); num_particles],
      particle_compliance: vec![placeholder_compliance(); num_particles],
      old,
    }
  }

  pub fn add_velocity_results(&mut self, results: HashMap<usize, VelocityTaskParticleData>) {
    for (id, data) in results {
      self.half_velocity[id] = data.half_velocity;
      self.particle_compliance[id] = data.compliance;

      let particle = &self.old.particles()[id];
      let mut new_particle = particle.reset_clone();
      new_particle.update_position(data.new_position);
      new_particle.set_thermostat_work(data.thermostat_work);
      self.local_container.add_particle(new_particle.clone());
      self.read_container.add_particle(new_particle);
    }
  }

  /// Replaces converting to a separate dense representation: `is_complete` is a single O(n)
  /// scan over the containers' own presence flags, with no extra allocation or particle copy.
  pub fn build(self) -> Option<IntegrationCache> {
    if !self.local_container.is_complete() || !self.read_container.is_complete() {
      return None;
    }

    Some(IntegrationCache::new(
      self.local_container,
      Arc::new(self.read_container),
      self.half_velocity,
      self.particle_compliance,
    ))
  }
}
