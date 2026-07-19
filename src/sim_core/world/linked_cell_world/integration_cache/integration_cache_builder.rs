use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::ParticleCompliance;
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::linked_cell_world::LinkedCellContainerOld;
use crate::sim_core::world::linked_cell_world::integration_cache::IntegrationCache;

pub struct IntegrationCacheBuilder {
  old: Arc<LinkedCellContainerOld>,
  local_container: LinkedCellContainerOld,
  read_container: LinkedCellContainerOld,
  half_velocity: Vec<Option<Vector3<f64>>>,
  particle_compliance: Vec<Option<ParticleCompliance>>,
}

impl IntegrationCacheBuilder {
  pub fn new(old: Arc<LinkedCellContainerOld>) -> Self {
    let num_particles = old.particles().len();
    let config = old.config().clone();
    let edge_condition = old.edge_condition();
    IntegrationCacheBuilder {
      local_container: LinkedCellContainerOld::new_empty(num_particles, config.clone(), edge_condition),
      read_container: LinkedCellContainerOld::new_empty(num_particles, config, edge_condition),
      half_velocity: vec![None; num_particles],
      particle_compliance: vec![None; num_particles],
      old,
    }
  }

  pub fn add_velocity_results(&mut self, results: HashMap<usize, VelocityTaskParticleData>) {
    for (id, data) in results {
      self.half_velocity[id] = Some(data.half_velocity);
      self.particle_compliance[id] = Some(data.compliance);

      let particle = self.old.particles()[id].as_ref().unwrap();
      let mut new_particle = particle.reset_clone();
      new_particle.update_position(data.new_position);
      new_particle.set_thermostat_work(data.thermostat_work);
      let arc = Arc::new(new_particle);
      self.local_container.add_particle(Arc::clone(&arc));
      self.read_container.add_particle(arc);
    }
  }

  pub fn build(self) -> Option<IntegrationCache> {
    if !self.half_velocity.iter().all(|o| o.is_some()) {
      return None;
    }

    let half_velocity_cache = self.half_velocity.into_iter().map(|o| o.unwrap()).collect();
    let particle_compliance = self.particle_compliance.into_iter().map(|o| o.unwrap()).collect();

    Some(IntegrationCache::new(
      self.local_container,
      Arc::new(self.read_container),
      half_velocity_cache,
      particle_compliance,
    ))
  }
}
