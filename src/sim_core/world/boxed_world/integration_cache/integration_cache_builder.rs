use crate::data::SimulationConfig;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::periodic::apply_velocity_constraint_periodic;
use crate::sim_core::world::boundary_constraint::{EdgeCondition, ParticleCompliance};
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;
use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;

pub struct IntegrationCacheBuilder {
  config: SimulationConfig,
  particles: HashMap<usize, Particle>,
  local_boxes: BoxContainer<Arc<SimulationBox>>,
  half_velocity: HashMap<usize, Vector3<f64>>,
  particle_compliance: HashMap<usize, ParticleCompliance>,
}

impl IntegrationCacheBuilder {
  pub fn new(
    config: SimulationConfig,
    box_container_config: BoxContainerConfig,
    particles: HashMap<usize, Particle>,
  ) -> Self {
    let num_particles = particles.len();
    IntegrationCacheBuilder {
      config,
      particles,
      local_boxes: BoxContainer::<Arc<SimulationBox>>::from_config(box_container_config),
      half_velocity: HashMap::with_capacity(num_particles),
      particle_compliance: HashMap::with_capacity(num_particles),
    }
  }

  pub fn add_velocity_results(&mut self, results: HashMap<usize, VelocityTaskParticleData>) {
    for (id, data) in results {
      let particle = self.particles.get_mut(&id).unwrap();
      particle.update_position(data.new_position);
      particle.set_thermostat_work(data.thermostat_work);

      self.local_boxes.add_particle(Arc::new(particle.clone()));

      self.half_velocity.insert(id, data.half_velocity);
      self.particle_compliance.insert(id, data.compliance);
    }
  }

  /// Builds IntegrationCache. Returns None if not all particles have received velocity results.
  pub fn build(self) -> Option<IntegrationCache> {
    let all_supplied = self.particles.len() == self.half_velocity.len()
      && self
        .particles
        .keys()
        .all(|id| self.half_velocity.contains_key(id));

    if !all_supplied {
      return None;
    }

    #[cfg(debug_assertions)]
    for sim_box in self.local_boxes.simulation_boxes().iter() {
      for particle in sim_box.as_ref().particles().values() {
        log::debug!(
          "Particle {} is stored in simulation box {}",
          particle.get_id(),
          sim_box.id()
        );
      }
    }

    Some(IntegrationCache::new(
      self.local_boxes,
      self.half_velocity,
      self.particle_compliance,
    ))
  }
}
