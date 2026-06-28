use std::sync::Arc;

use crate::data::types::AtomType;
use crate::data::{ParticleConfig, SimulationConfig};
use crate::particle::Particle;
use crate::persistence::dto::world::boxed::BoxedWorldDTO;
use crate::persistence::dto::world::boxed::BoxedWorldDTOWithoutHistory;
use crate::persistence::dto::world::history::HistoryDTO;
use crate::sim_core::world::boxed_world::box_container::box_container_config;
use crate::sim_core::world::linked_cell_world::LinkedCellContainer;

mod getter;

pub struct LinkedCellHistoryManager {
  config: SimulationConfig,
  thermostat_epsilon: Vec<f64>,
  history: Vec<Arc<LinkedCellContainer>>,
  current_index: usize,
}

impl LinkedCellHistoryManager {
  pub fn with_config(config: SimulationConfig, particle_config: ParticleConfig) -> Self {
    let atoms: Vec<Particle> = particle_config.atoms;
    let container_config = box_container_config::new_config(&atoms, config.world_size);
    let particles: Vec<Arc<Particle>> = atoms.into_iter().map(Arc::new).collect();

    let mut container = LinkedCellContainer::new(particles, container_config, config.edge_condition);
    container.sort();

    let mut thermostat_epsilon = Vec::with_capacity(config.max_iteration_till_reset);
    thermostat_epsilon.push(0.);

    let mut history = Vec::with_capacity(config.max_iteration_till_reset);
    history.push(Arc::new(container));

    LinkedCellHistoryManager {
      config,
      thermostat_epsilon,
      history,
      current_index: 0,
    }
  }

  pub fn current_container(&self) -> Arc<LinkedCellContainer> {
    self.history.last().unwrap().clone()
  }

  pub fn current_index(&self) -> usize {
    self.current_index
  }

  pub fn reset_clone(&self) -> LinkedCellHistoryManager {
    let mut thermostat_epsilon = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    thermostat_epsilon.push(*self.thermostat_epsilon.last().unwrap());

    let mut history = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    history.push(self.history.last().unwrap().clone());

    LinkedCellHistoryManager {
      config: self.config.clone(),
      thermostat_epsilon,
      history,
      current_index: 0,
    }
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

  pub fn push_container(&mut self, container: LinkedCellContainer) {
    self.history.push(Arc::new(container));
    self.current_index += 1;
  }

  pub fn reset_container(&mut self) {
    let mut new_thermostat_epsilon = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    if let Some(last_epsilon) = self.thermostat_epsilon.pop() {
      new_thermostat_epsilon.push(last_epsilon);
    } else {
      panic!("Thermostat epsilon is empty!");
    }

    let mut new_history = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    if let Some(last_container) = self.history.pop() {
      new_history.push(last_container);
    } else {
      panic!("History is empty!");
    }

    self.current_index = 0;
    self.history = new_history;
    self.thermostat_epsilon = new_thermostat_epsilon;
  }

  pub fn to_dto(self, partial: BoxedWorldDTOWithoutHistory, lower_index: usize) -> BoxedWorldDTO {
    let history = self.to_history_dto(lower_index);
    partial.with_history(history)
  }

  fn to_history_dto(self, lower_index: usize) -> HistoryDTO {
    let box_container = self.history[lower_index..]
      .iter()
      .map(|container| container.to_transfer_struct())
      .collect();
    HistoryDTO { box_container, thermostat_epsilon: self.thermostat_epsilon }
  }

  pub fn to_transfer_struct(&self, lower_index: usize) -> HistoryDTO {
    let box_container = self.history[lower_index..]
      .iter()
      .map(|container| container.to_transfer_struct())
      .collect();
    HistoryDTO { box_container, thermostat_epsilon: self.thermostat_epsilon.clone() }
  }

  pub fn get_particle_counts(&self) -> (usize, usize, usize) {
    let mut c_count = 0;
    let mut fe_count = 0;

    if let Some(container) = self.history.last() {
      for particle in container.particles().iter().filter_map(|p| p.as_ref()) {
        match particle.get_type() {
          AtomType::C | AtomType::C_nanotube => c_count += 1,
          AtomType::Fe => fe_count += 1,
        }
      }
    }

    let total = c_count + fe_count;
    (total, c_count, fe_count)
  }
}
