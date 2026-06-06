use crate::data::SimulationConfig;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::ParticleCompliance;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;
use nalgebra::Vector3;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

pub struct IntegrationCacheBuilder {
  config: SimulationConfig,
  particles: HashMap<usize, Particle>,
  local_boxes: BoxContainer<SimulationBox>,
  half_velocity: HashMap<usize, Vector3<f64>>,
  particle_compliance: HashMap<usize, ParticleCompliance>,
  pair_corrected_ids: HashSet<usize>,
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
      local_boxes: BoxContainer::<SimulationBox>::new_local(box_container_config),
      half_velocity: HashMap::with_capacity(num_particles),
      particle_compliance: HashMap::with_capacity(num_particles),
      pair_corrected_ids: HashSet::new(),
    }
  }

  pub fn add_velocity_results(
    &mut self,
    results: HashMap<usize, VelocityTaskParticleData>,
    pair_corrected: bool,
  ) {
    if pair_corrected {
      for id in results.keys() {
        self.pair_corrected_ids.insert(*id);
      }
    }

    for (id, data) in results {
      let mut particle = self.particles.get(&id).unwrap().clone();
      particle.update_position(data.new_position);
      particle.set_thermostat_work(data.thermostat_work);

      self.local_boxes.add_particle(Arc::new(particle));

      self.half_velocity.insert(id, data.half_velocity);
      self.particle_compliance.insert(id, data.compliance);
    }
  }

  /// Builds IntegrationCache. Returns None if not all particles have received velocity results.
  pub fn build(self) -> Option<IntegrationCache> {
    let all_supplied = !self.particles.is_empty()
      && self.particles.len() == self.half_velocity.len()
      && self
        .particles
        .keys()
        .all(|id| self.half_velocity.contains_key(id));

    if !all_supplied {
      return None;
    }

    let shared_boxes = self.local_boxes.into_shared();

    #[cfg(debug_assertions)]
    for sim_box in shared_boxes.simulation_boxes().iter() {
      for particle in sim_box.particles().values() {
        log::debug!(
          "Particle {} is stored in simulation box {}",
          particle.get_id(),
          sim_box.id()
        );
      }
    }

    Some(IntegrationCache::new(
      shared_boxes,
      self.half_velocity,
      self.particle_compliance,
      self.pair_corrected_ids,
    ))
  }
}
