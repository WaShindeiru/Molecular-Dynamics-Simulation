use crate::data::types::AtomType;
use crate::data::{ParticleConfig, SimulationConfig};
use crate::persistence::dto::world::history::HistoryDTO;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use std::sync::Arc;

mod getter;

pub struct HistoryManager {
  config: SimulationConfig,
  thermostat_epsilon: Vec<f64>,
  box_container_config: BoxContainerConfig,
  history: Vec<Arc<BoxContainer<Arc<SimulationBox>>>>,
  current_index: usize,
}
impl HistoryManager {
  pub fn with_config(config: SimulationConfig, particle_config: ParticleConfig) -> Self {
    let atoms_to_use: Vec<Particle> = particle_config.atoms;

    let mut thermostat_epsilon: Vec<f64> = Vec::with_capacity(config.max_iteration_till_reset);
    thermostat_epsilon.push(0.);

    let mut history: Vec<Arc<BoxContainer>> = Vec::with_capacity(config.max_iteration_till_reset);
    history.push(Arc::new(BoxContainer::new(atoms_to_use, config.world_size)));

    let box_container_config = history.get(0).unwrap().config().clone();

    HistoryManager {
      config,
      thermostat_epsilon,
      box_container_config,
      history,
      current_index: 0,
    }
  }

  pub fn current_box_container(&self) -> Arc<BoxContainer<Arc<SimulationBox>>> {
    self.history.last().unwrap().clone()
  }

  pub fn current_index(&self) -> usize {
    self.current_index
  }

  pub fn clone_for_cache(&self) -> (Vec<Arc<BoxContainer<Arc<SimulationBox>>>>, Vec<f64>) {
    (self.history.clone(), self.thermostat_epsilon.clone())
  }

  pub fn current_thermostat_epsilon(&self) -> f64 {
    *self.thermostat_epsilon.last().unwrap()
  }

  pub fn thermostat_epsilon_of_iteration(&self, iteration: usize) -> f64 {
    *self.thermostat_epsilon.get(iteration).unwrap()
  }

  pub fn add_thermostat_epsilon(&mut self, thermostat_epsilon: f64) {
    self.thermostat_epsilon.push(thermostat_epsilon);
  }

  // TODO: should I update current_index here?
  pub fn push_box_container(&mut self, container: BoxContainer<Arc<SimulationBox>>) {
    self.history.push(Arc::new(container));
    self.current_index += 1;
  }

  pub fn reset_container(&mut self) {
    let new_index = 0;

    let mut new_thermostat_epsilon: Vec<f64> =
      Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    if let Some(last_epsilon) = self.thermostat_epsilon.pop() {
      new_thermostat_epsilon.push(last_epsilon);
    } else {
      panic!("Thermostat epsilon is empty!");
    }

    let mut new_box_container: Vec<Arc<BoxContainer>> =
      Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    if let Some(last_history) = self.history.pop() {
      new_box_container.push(last_history);
    } else {
      panic!("Boxes are empty!");
    }

    self.current_index = new_index;
    self.history = new_box_container;
    self.thermostat_epsilon = new_thermostat_epsilon;
  }

  pub fn to_transfer_struct(&self, lower_index: usize) -> HistoryDTO {
    let box_container = self.history[lower_index..]
      .iter()
      .map(|container| container.to_transfer_struct())
      .collect();

    HistoryDTO {
      box_container,
      thermostat_epsilon: self.thermostat_epsilon.clone(),
    }
  }

  pub fn get_particle_counts(&self) -> (usize, usize, usize) {
    let mut c_count = 0;
    let mut fe_count = 0;

    if let Some(container) = self.history.last() {
      for sim_box in container.simulation_boxes().iter() {
        for (_, particle) in sim_box.particles() {
          match particle.get_type() {
            AtomType::C => c_count += 1,
            AtomType::Fe => fe_count += 1,
          }
        }
      }
    }

    let total = c_count + fe_count;
    (total, c_count, fe_count)
  }
}
