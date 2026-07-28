mod getter;

use std::sync::Arc;

use crate::data::types::AtomType;
use crate::data::{ParticleConfig, SimulationConfig};
use crate::particle::Particle;
use crate::persistence::dto::world::boxed::BoxedWorldDTO;
use crate::persistence::dto::world::boxed::BoxedWorldDTOWithoutHistory;
use crate::persistence::dto::world::history::HistoryDTO;
use crate::sim_core::world::cell::{LinkedCellContainer, box_container_config};
use crate::sim_core::world::thermostat::IntegrationAlgorithm;

fn has_nanotube_thermostat(config: &SimulationConfig) -> bool {
  matches!(
    config.integration_algorithm,
    IntegrationAlgorithm::NoseHooverVerlet { nanotube_thermostat: Some(_), .. }
  )
}

pub struct OptimizedHistoryManager {
  config: SimulationConfig,
  thermostat_epsilon: Vec<f64>,
  temperature: Vec<f64>,
  nanotube_thermostat_epsilon: Option<Vec<f64>>,
  nanotube_temperature: Option<Vec<f64>>,
  history: Vec<Arc<LinkedCellContainer>>,
  current_index: usize,
}

impl OptimizedHistoryManager {
  pub fn with_config(config: SimulationConfig, particle_config: ParticleConfig) -> Self {
    let atoms: Vec<Particle> = particle_config.atoms;
    let container_config = box_container_config::new_config(&atoms, config.world_size);

    let mut container = LinkedCellContainer::new(atoms, container_config, config.edge_condition);
    container.sort();

    let mut thermostat_epsilon = Vec::with_capacity(config.max_iteration_till_reset);
    thermostat_epsilon.push(0.);

    let mut temperature = Vec::with_capacity(config.max_iteration_till_reset);
    temperature.push(0.);

    let nanotube_active = has_nanotube_thermostat(&config);
    let nanotube_thermostat_epsilon = nanotube_active.then(|| vec![0.]);
    let nanotube_temperature = nanotube_active.then(|| vec![0.]);

    let mut history = Vec::with_capacity(config.max_iteration_till_reset);
    history.push(Arc::new(container));

    OptimizedHistoryManager {
      config,
      thermostat_epsilon,
      temperature,
      nanotube_thermostat_epsilon,
      nanotube_temperature,
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

  pub fn reset_clone(&self) -> OptimizedHistoryManager {
    let mut thermostat_epsilon = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    thermostat_epsilon.push(*self.thermostat_epsilon.last().unwrap());

    let mut temperature = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    temperature.push(*self.temperature.last().unwrap());

    let carry_last = |v: &Option<Vec<f64>>| v.as_ref().map(|vec| vec![*vec.last().unwrap()]);
    let nanotube_thermostat_epsilon = carry_last(&self.nanotube_thermostat_epsilon);
    let nanotube_temperature = carry_last(&self.nanotube_temperature);

    let mut history = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    history.push(self.history.last().unwrap().clone());

    OptimizedHistoryManager {
      config: self.config.clone(),
      thermostat_epsilon,
      temperature,
      nanotube_thermostat_epsilon,
      nanotube_temperature,
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

  pub fn add_temperature(&mut self, temperature: f64) {
    self.temperature.push(temperature);
  }

  /// `None` unless the nanotube thermostat is configured, in which case `epsilon` must be `Some`.
  pub fn current_nanotube_thermostat_epsilon(&self) -> Option<f64> {
    self.nanotube_thermostat_epsilon.as_ref().map(|v| *v.last().unwrap())
  }

  pub fn add_nanotube_thermostat_epsilon(&mut self, epsilon: Option<f64>) {
    if let Some(vec) = &mut self.nanotube_thermostat_epsilon {
      vec.push(epsilon.expect("nanotube thermostat is configured but no epsilon was computed"));
    }
  }

  pub fn add_nanotube_temperature(&mut self, temperature: Option<f64>) {
    if let Some(vec) = &mut self.nanotube_temperature {
      vec.push(temperature.expect("nanotube thermostat is configured but no temperature was computed"));
    }
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

    let mut new_temperature = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    if let Some(last_temperature) = self.temperature.pop() {
      new_temperature.push(last_temperature);
    } else {
      panic!("Temperature is empty!");
    }

    let pop_last = |v: &mut Option<Vec<f64>>| -> Option<Vec<f64>> {
      v.as_mut().map(|vec| vec![vec.pop().expect("nanotube history vec is empty")])
    };
    let new_nanotube_thermostat_epsilon = pop_last(&mut self.nanotube_thermostat_epsilon);
    let new_nanotube_temperature = pop_last(&mut self.nanotube_temperature);

    let mut new_history = Vec::with_capacity(self.config.max_iteration_till_reset + 1);
    if let Some(last_container) = self.history.pop() {
      new_history.push(last_container);
    } else {
      panic!("History is empty!");
    }

    self.current_index = 0;
    self.history = new_history;
    self.thermostat_epsilon = new_thermostat_epsilon;
    self.temperature = new_temperature;
    self.nanotube_thermostat_epsilon = new_nanotube_thermostat_epsilon;
    self.nanotube_temperature = new_nanotube_temperature;
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
    HistoryDTO {
      box_container,
      thermostat_epsilon: self.thermostat_epsilon,
      temperature: self.temperature,
      nanotube_thermostat_epsilon: self.nanotube_thermostat_epsilon,
      nanotube_temperature: self.nanotube_temperature,
    }
  }

  pub fn to_transfer_struct(&self, lower_index: usize) -> HistoryDTO {
    let box_container = self.history[lower_index..]
      .iter()
      .map(|container| container.to_transfer_struct())
      .collect();
    HistoryDTO {
      box_container,
      thermostat_epsilon: self.thermostat_epsilon.clone(),
      temperature: self.temperature.clone(),
      nanotube_thermostat_epsilon: self.nanotube_thermostat_epsilon.clone(),
      nanotube_temperature: self.nanotube_temperature.clone(),
    }
  }

  pub fn get_particle_counts(&self) -> (usize, usize, usize) {
    let mut c_count = 0;
    let mut fe_count = 0;

    if let Some(container) = self.history.last() {
      for particle in container.particles().iter() {
        match particle.get_type() {
          AtomType::C | AtomType::C_nanotube | AtomType::C_nanotube_static => c_count += 1,
          AtomType::Fe => fe_count += 1,
        }
      }
    }

    let total = c_count + fe_count;
    (total, c_count, fe_count)
  }
}
