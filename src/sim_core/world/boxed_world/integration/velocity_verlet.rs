use std::sync::Arc;

use crate::sim_core::world::boxed_world::BoxedWorld;
use crate::sim_core::world::thermostat::IntegrationStateUpdateResponse;

impl BoxedWorld {
  pub fn update_velocity_verlet(&mut self, next_iteration: usize) {
    let current_thermostat_epsilon = self
      .persistance_reset
      .history_manager()
      .current_thermostat_epsilon();
    let current_box_container = self
      .persistance_reset
      .history_manager()
      .current_box_container();

    let integration_cache = self.task_manager.half_velocity_step(
      Arc::clone(&current_box_container),
      current_thermostat_epsilon,
      next_iteration,
    );

    self.integration_cache = Some(Arc::clone(&integration_cache));

    let new_thermostat_epsilon = 0.;

    self
      .persistance_reset
      .history_manager_mut()
      .add_thermostat_epsilon(new_thermostat_epsilon);

    let mut computation_collector = self.task_manager.force_step(Arc::clone(&integration_cache));

    computation_collector.apply_gravity(next_iteration);

    computation_collector.set_velocity(new_thermostat_epsilon);

    let current_custom_velocities =
      self.velocity_manager.compute_velocities_for_iteration(next_iteration);
    computation_collector.apply_custom_velocities(&current_custom_velocities);

    let simulation_temperature = computation_collector.get_mean_temperature();

    let result = self.integration_algorithm_state.update_state(
      next_iteration,
      self.config.time_step,
      &self.config.integration_algorithm,
      simulation_temperature,
      computation_collector.particles(),
    );
    match result {
      IntegrationStateUpdateResponse::VelocityVerlet => {}
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
