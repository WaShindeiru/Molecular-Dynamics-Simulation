use log::info;

use crate::data::units::TEMPERATURE_U;
use crate::sim_core::world::optimized_world::OptimizedWorld;
use crate::sim_core::world::optimized_world::integration::verlet_nose_hoover::thermostat::compute_new_thermostat_epsilon;
use crate::sim_core::world::thermostat::{
  IntegrationAlgorithm, IntegrationAlgorithmState, IntegrationStateUpdateResponse,
};

impl OptimizedWorld {
  pub fn update_verlet_nose_hoover(&mut self, next_iteration: usize) {
    let current_thermostat_epsilon = self
      .persistance_reset
      .history_manager()
      .current_thermostat_epsilon();
    let current_container = self
      .persistance_reset
      .history_manager()
      .current_container();

    let temperature_index = if let IntegrationAlgorithmState::NoseHooverVerlet {
      temperature_index, ..
    } = self.integration_algorithm_state
    {
      temperature_index
    } else {
      panic!("Expected NoseHooverVerlet integration state!")
    };

    let (current_desired_temperature, q_effective_mass) =
      if let IntegrationAlgorithm::NoseHooverVerlet { desired_temperature, q_effective_mass: q } =
        &self.config.integration_algorithm
      {
        let temp_info = desired_temperature.get(temperature_index).unwrap();
        (temp_info.desired_temperature, *q)
      } else {
        panic!("Expected NoseHooverVerlet integration!")
      };

    let integration_cache = self.task_manager.half_velocity_step(
      current_container,
      current_thermostat_epsilon,
      next_iteration,
      self.config.time_step,
    );

    let new_thermostat_epsilon = compute_new_thermostat_epsilon(
      current_thermostat_epsilon,
      &integration_cache.half_velocity_cache,
      integration_cache
        .local_container
        .particles()
        .iter()
        .filter(|p| p.is_atom()),
      self.config.time_step,
      q_effective_mass,
      current_desired_temperature,
    );

    self
      .persistance_reset
      .history_manager_mut()
      .add_thermostat_epsilon(new_thermostat_epsilon);

    let mut computation_collector = self.task_manager.force_step(integration_cache);

    computation_collector.apply_gravity(next_iteration);
    computation_collector.set_velocity(new_thermostat_epsilon);

    let current_custom_velocities =
      self.velocity_manager.compute_velocities_for_iteration(next_iteration);
    computation_collector.apply_custom_velocities(&current_custom_velocities);

    let current_control_velocities = self
      .control_velocity_manager
      .compute_controlled_velocities_for_iteration(next_iteration);
    computation_collector.compute_controlled_velocity_particles(&current_control_velocities, self.config.alpha);

    let simulation_temperature = computation_collector.get_mean_temperature();

    let result = self.integration_algorithm_state.update_state(
      next_iteration,
      self.config.time_step,
      &self.config.integration_algorithm,
      simulation_temperature,
      computation_collector.particles(),
    );

    match result {
      IntegrationStateUpdateResponse::NoseHooverVerlet {
        updated,
        temperature,
      } => {
        if updated {
          let entry = self
            .integration_algorithm_state
            .get_previous_history_entry()
            .unwrap();

          let temperature_started = entry.temperature_started.unwrap().temperature;
          let temperature_achieved = entry.temperature_achieved.unwrap().temperature;
          let temperature_switched = entry.temperature_switched.unwrap().temperature;

          info!(
            "Iteration: {iteration}, current simulation temperature: {simulation_temperature} K, \
             temperature {achieved_temperature} K achieved, switching to {next_temperature} K.\
             temperature_started: {temperature_started} K, temperature_achieved: {temperature_achieved} K, \
             temperature_switched: {temperature_switched} K",
            iteration = self.iteration,
            simulation_temperature = simulation_temperature * TEMPERATURE_U,
            achieved_temperature = current_desired_temperature * TEMPERATURE_U,
            next_temperature = temperature * TEMPERATURE_U,
            temperature_started = temperature_started * TEMPERATURE_U,
            temperature_achieved = temperature_achieved * TEMPERATURE_U,
            temperature_switched = temperature_switched * TEMPERATURE_U,
          );
        }
      }
      _ => panic!("Wrong result type"),
    }

    self
      .persistance_reset
      .history_manager_mut()
      .push_container(computation_collector.build());

    self.iteration += 1;
    assert_eq!(self.iteration, next_iteration);
  }
}
