use std::collections::{HashMap, HashSet};
use std::sync::Arc;

use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;


pub struct ComputationCollector {
  particles_modified: HashMap<usize, Particle>,
  particle_box_mapping: HashMap<usize, usize>,
  box_container_config: BoxContainerConfig,
  force_visited: HashSet<usize>,
}

impl ComputationCollector {
  pub fn from_integration_cache(cache: &IntegrationCache) -> Self {
    let particles_modified = cache.integration_box_cache.all_particles_cloned();
    let particle_box_mapping = cache.integration_box_cache.box_id_cache().clone();
    let box_container_config = *cache.integration_box_cache.config();

    ComputationCollector {
      particles_modified,
      particle_box_mapping,
      box_container_config,
      force_visited: HashSet::new(),
    }
  }

  pub fn into_box_container(self) -> BoxContainer<Arc<SimulationBox>> {
    BoxContainer::from_particles(
      self.box_container_config,
      &self.particles_modified,
      &self.particle_box_mapping,
    )
  }
}