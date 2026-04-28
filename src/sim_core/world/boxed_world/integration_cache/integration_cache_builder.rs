use std::collections::HashMap;
use std::sync::Arc;
use nalgebra::Vector3;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::ParticleCompliance;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;

pub struct IntegrationCacheBuilder {
  particles: HashMap<usize, Particle>,
  local_boxes: BoxContainer<SimulationBox>,
  half_velocity: HashMap<usize, Vector3<f64>>,
  particle_compliance: HashMap<usize, ParticleCompliance>,
}

impl IntegrationCacheBuilder {
  pub fn new(box_container_config: BoxContainerConfig, particles: HashMap<usize, Particle>) -> Self {
    let num_particles = particles.len();
    IntegrationCacheBuilder {
      particles,
      local_boxes: BoxContainer::<SimulationBox>::new_local(box_container_config),
      half_velocity: HashMap::with_capacity(num_particles),
      particle_compliance: HashMap::with_capacity(num_particles),
    }
  }

  pub fn add_velocity_results(&mut self, results: HashMap<usize, VelocityTaskParticleData>) {
    for (id, data) in results {
      let particle = self.particles.get_mut(&id).unwrap();
      particle.update_position(data.new_position);
      particle.set_thermostat_work(data.thermostat_work);

      let box_id = self.local_boxes.config().box_id_for_position(particle.get_position());
      self.local_boxes.get_box_mut(box_id).add_particle(Arc::new(particle.clone()));

      self.half_velocity.insert(id, data.half_velocity);
      self.particle_compliance.insert(id, data.compliance);
    }
  }

  /// Builds IntegrationCache. Returns None if not all particles have received velocity results.
  pub fn build(self) -> Option<IntegrationCache> {
    let all_supplied = self.particles.len() == self.half_velocity.len()
      && self.particles.keys().all(|id| self.half_velocity.contains_key(id));

    if !all_supplied {
      return None;
    }

    Some(IntegrationCache::new(
      self.local_boxes.into_shared(),
      self.half_velocity,
      self.particle_compliance,
    ))
  }
}
