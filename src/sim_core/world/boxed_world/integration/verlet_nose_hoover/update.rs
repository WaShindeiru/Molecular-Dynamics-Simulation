use std::sync::Arc;

use log::info;

use crate::data::units::TEMPERATURE_U;

use crate::sim_core::world::boxed_world::BoxedWorld;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::thermostat::compute_new_thermostat_epsilon;

use crate::sim_core::world::integration::IntegrationAlgorithm;
use crate::sim_core::world::integration::IntegrationAlgorithmState;
use crate::sim_core::world::integration::IntegrationStateUpdateResponse;

impl BoxedWorld {
  pub fn update_verlet_nose_hoover(&mut self, next_iteration: usize) {
    let current_thermostat_epsilon = self
      .persistance_reset
      .history_manager()
      .current_thermostat_epsilon();
    let current_box_container = self
      .persistance_reset
      .history_manager()
      .current_box_container();

    let current_desired_temperature: f64;
    let q_effective_mass: f64;

    let temperature_index = if let IntegrationAlgorithmState::NoseHooverVerlet {
      temperature_index,
      ..
    } = self.integration_algorithm_state
    {
      temperature_index
    } else {
      panic!("Expected NoseHooverVerlet integration state!")
    };

    if let IntegrationAlgorithm::NoseHooverVerlet {
      desired_temperature,
      q_effective_mass: q,
    } = &self.config.integration_algorithm
    {
      let temp_info = desired_temperature.get(temperature_index).unwrap();
      current_desired_temperature = temp_info.desired_temperature;
      q_effective_mass = *q;
    } else {
      panic!("Expected NoseHooverVerlet integration!")
    }

    let integration_cache = self.task_manager.half_velocity_step(
      Arc::clone(&current_box_container),
      current_thermostat_epsilon,
      next_iteration,
    );

    self.integration_cache = Some(Arc::clone(&integration_cache));

    let new_thermostat_epsilon;
    {
      new_thermostat_epsilon = compute_new_thermostat_epsilon(
        current_thermostat_epsilon,
        integration_cache.half_velocity_cache(),
        integration_cache.box_cache().all_particles(),
        self.config.time_step,
        q_effective_mass,
        current_desired_temperature,
      )
    }

    self
      .persistance_reset
      .history_manager_mut()
      .add_thermostat_epsilon(new_thermostat_epsilon);

    let mut computation_collector = self.task_manager.force_step(Arc::clone(&integration_cache));

    computation_collector.apply_gravity();

    computation_collector.set_velocity(new_thermostat_epsilon);

    let simulation_temperature = computation_collector.get_mean_temperature();

    let result = self.integration_algorithm_state.update_state(
      next_iteration,
      self.config.time_step,
      &self.config.integration_algorithm,
      simulation_temperature,
    );
    match result {
      IntegrationStateUpdateResponse::NoseHooverVerlet {
        updated,
        temperature,
      } => {
        if updated {
          info!(
            "Iteration: {}, Simulation temperature: {}, Temperature: {} achieved, switching to temperature: {}",
            self.iteration,
            simulation_temperature * TEMPERATURE_U,
            current_desired_temperature * TEMPERATURE_U,
            temperature * TEMPERATURE_U
          );
        }
      }
      _ => panic!("Wrong result type"),
    }

    self
      .persistance_reset
      .history_manager_mut()
      .push_box_container(computation_collector.into_box_container());

    self.iteration += 1;
    assert_eq!(self.iteration, next_iteration);
  }
}
