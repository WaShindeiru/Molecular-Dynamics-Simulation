use log::info;
use crate::data::units::TEMPERATURE_U;
use crate::sim_core::world::optimized_world::OptimizedWorld;
use crate::sim_core::world::optimized_world::computation_collector::ComputationCollector;
use crate::sim_core::world::optimized_world::integration::verlet_nose_hoover::thermostat::compute_new_thermostat_epsilon;
use crate::sim_core::world::thermostat::{
  IntegrationAlgorithm, IntegrationAlgorithmState, IntegrationStateUpdateResponse,
};

impl OptimizedWorld {
  fn current_desired_temperature_and_q(&self) -> (f64, f64) {
    let temperature_index = if let IntegrationAlgorithmState::NoseHooverVerlet {
      temperature_index, ..
    } = self.integration_algorithm_state
    {
      temperature_index
    } else {
      panic!("Expected NoseHooverVerlet integration state!")
    };

    if let IntegrationAlgorithm::NoseHooverVerlet { desired_temperature, q_effective_mass: q, .. } =
      &self.config.integration_algorithm
    {
      let temp_info = desired_temperature.get(temperature_index).unwrap();
      (temp_info.desired_temperature, *q)
    } else {
      panic!("Expected NoseHooverVerlet integration!")
    }
  }

  /// `None` unless a nanotube thermostat is configured.
  fn current_nanotube_desired_temperature_and_q(&self) -> Option<(f64, f64)> {
    let Some(thermostat_state) = &self.thermostat_integration_algorithm_state else {
      return None;
    };

    let temperature_index = if let IntegrationAlgorithmState::NoseHooverVerlet {
      temperature_index, ..
    } = thermostat_state
    {
      *temperature_index
    } else {
      panic!("Expected NoseHooverVerlet integration state for nanotube thermostat!")
    };

    let IntegrationAlgorithm::NoseHooverVerlet { nanotube_thermostat: Some(nanotube_thermostat), .. } =
      &self.config.integration_algorithm
    else {
      panic!("thermostat_integration_algorithm_state is set but nanotube_thermostat is not configured!")
    };

    let temp_info = nanotube_thermostat.desired_temperature.get(temperature_index).unwrap();
    Some((temp_info.desired_temperature, nanotube_thermostat.q_effective_mass))
  }

  pub fn update_verlet_nose_hoover(&mut self, next_iteration: usize) {
    let current_thermostat_epsilon = self
      .persistance_reset
      .history_manager()
      .current_thermostat_epsilon();
    let current_nanotube_thermostat_epsilon = self
      .persistance_reset
      .history_manager()
      .current_nanotube_thermostat_epsilon();
    let current_container = self
      .persistance_reset
      .history_manager()
      .current_container();

    let (current_desired_temperature, q_effective_mass) = self.current_desired_temperature_and_q();
    let nanotube_desired_temperature_and_q = self.current_nanotube_desired_temperature_and_q();

    let integration_cache = self.task_manager.half_velocity_step(
      current_container,
      current_thermostat_epsilon,
      current_nanotube_thermostat_epsilon,
      next_iteration,
      self.config.time_step,
    );

    let (new_thermostat_epsilon, new_nanotube_thermostat_epsilon) =
      if let Some((nanotube_desired_temperature, nanotube_q_effective_mass)) =
        nanotube_desired_temperature_and_q
      {
        let new_thermostat_epsilon = compute_new_thermostat_epsilon(
          current_thermostat_epsilon,
          &integration_cache.half_velocity_cache,
          integration_cache
            .local_container
            .particles()
            .iter()
            .filter(|p| p.is_atom() && !p.is_nanotube_atom()),
          self.config.time_step,
          q_effective_mass,
          current_desired_temperature,
        );

        let new_nanotube_thermostat_epsilon = compute_new_thermostat_epsilon(
          current_nanotube_thermostat_epsilon.unwrap(),
          &integration_cache.half_velocity_cache,
          integration_cache.local_container.particles().iter().filter(|p| p.is_nanotube_atom()),
          self.config.time_step,
          nanotube_q_effective_mass,
          nanotube_desired_temperature,
        );

        (new_thermostat_epsilon, Some(new_nanotube_thermostat_epsilon))
      } else {
        let new_thermostat_epsilon = compute_new_thermostat_epsilon(
          current_thermostat_epsilon,
          &integration_cache.half_velocity_cache,
          integration_cache.local_container.particles().iter().filter(|p| p.is_atom()),
          self.config.time_step,
          q_effective_mass,
          current_desired_temperature,
        );

        (new_thermostat_epsilon, None)
      };

    self
      .persistance_reset
      .history_manager_mut()
      .add_thermostat_epsilon(new_thermostat_epsilon);
    self
      .persistance_reset
      .history_manager_mut()
      .add_nanotube_thermostat_epsilon(new_nanotube_thermostat_epsilon);

    let mut computation_collector = self.task_manager.force_step(integration_cache);

    computation_collector.apply_gravity(next_iteration);
    computation_collector.set_velocity(new_thermostat_epsilon, new_nanotube_thermostat_epsilon);

    let current_custom_velocities =
      self.velocity_manager.compute_velocities_for_iteration(next_iteration);
    computation_collector.apply_custom_velocities(&current_custom_velocities);

    let current_control_velocities = self
      .control_velocity_manager
      .compute_controlled_velocities_for_iteration(next_iteration);
    computation_collector.compute_controlled_velocity_particles(&current_control_velocities, self.config.alpha);

    self.update_integration_states(next_iteration, &computation_collector, current_desired_temperature);

    self
      .persistance_reset
      .history_manager_mut()
      .push_container(computation_collector.build());

    self.iteration += 1;
    assert_eq!(self.iteration, next_iteration);
  }

  /// Advances the main Nose-Hoover schedule and, if configured, the independent nanotube one -
  /// each drives its own [`IntegrationAlgorithmState`] and logs on its own switch.
  fn update_integration_states(
    &mut self,
    next_iteration: usize,
    computation_collector: &ComputationCollector,
    current_desired_temperature: f64,
  ) {
    let simulation_temperature = computation_collector.get_mean_temperature();

    let result = self.integration_algorithm_state.update_state(
      next_iteration,
      self.config.time_step,
      &self.config.integration_algorithm,
      simulation_temperature,
      computation_collector.particles(),
    );

    Self::log_if_updated(
      self.iteration,
      simulation_temperature,
      current_desired_temperature,
      result,
      &self.integration_algorithm_state,
      "[general] ",
    );

    self
      .persistance_reset
      .history_manager_mut()
      .add_temperature(simulation_temperature);

    if self.thermostat_integration_algorithm_state.is_none() {
      return;
    }

    // Read everything needed from `self` up front, as plain values, before taking the `&mut`
    // borrow of `thermostat_integration_algorithm_state` below - a method call needs all of
    // `&self`, which would otherwise conflict with that field's mutable borrow.
    let (nanotube_current_desired_temperature, _) = self
      .current_nanotube_desired_temperature_and_q()
      .expect("thermostat_integration_algorithm_state is set but nanotube_thermostat is not configured!");
    let IntegrationAlgorithm::NoseHooverVerlet { nanotube_thermostat: Some(nanotube_thermostat), .. } =
      &self.config.integration_algorithm
    else {
      panic!("thermostat_integration_algorithm_state is set but nanotube_thermostat is not configured!")
    };
    let nanotube_algorithm = nanotube_thermostat.as_integration_algorithm();
    let time_step = self.config.time_step;
    let iteration = self.iteration;

    let nanotube_simulation_temperature = computation_collector.get_nanotube_mean_temperature();

    let thermostat_state = self.thermostat_integration_algorithm_state.as_mut().unwrap();

    let nanotube_result = thermostat_state.update_state(
      next_iteration,
      time_step,
      &nanotube_algorithm,
      nanotube_simulation_temperature,
      computation_collector.particles().filter(|p| p.is_nanotube_atom()),
    );

    Self::log_if_updated(
      iteration,
      nanotube_simulation_temperature,
      nanotube_current_desired_temperature,
      nanotube_result,
      thermostat_state,
      "[nanotube] ",
    );

    self
      .persistance_reset
      .history_manager_mut()
      .add_nanotube_temperature(Some(nanotube_simulation_temperature));
  }

  fn log_if_updated(
    iteration: usize,
    simulation_temperature: f64,
    current_desired_temperature: f64,
    result: IntegrationStateUpdateResponse,
    state: &IntegrationAlgorithmState,
    log_prefix: &str,
  ) {
    let IntegrationStateUpdateResponse::NoseHooverVerlet { updated, temperature } = result else {
      panic!("Wrong result type")
    };

    if !updated {
      return;
    }

    // Only looked up once a switch has actually happened - safe here since
    // `temperature_index` has already advanced past 0 by the time `updated` is true.
    let entry = state.get_previous_history_entry().unwrap();

    let temperature_started = entry.temperature_started.unwrap().temperature;
    let temperature_achieved = entry.temperature_achieved.unwrap().temperature;
    let temperature_switched = entry.temperature_switched.unwrap().temperature;

    info!(
      "{log_prefix}Iteration: {iteration}, current simulation temperature: {simulation_temperature} K, \
       temperature {achieved_temperature} K achieved, switching to {next_temperature} K.\
       temperature_started: {temperature_started} K, temperature_achieved: {temperature_achieved} K, \
       temperature_switched: {temperature_switched} K",
      simulation_temperature = simulation_temperature * TEMPERATURE_U,
      achieved_temperature = current_desired_temperature * TEMPERATURE_U,
      next_temperature = temperature * TEMPERATURE_U,
      temperature_started = temperature_started * TEMPERATURE_U,
      temperature_achieved = temperature_achieved * TEMPERATURE_U,
      temperature_switched = temperature_switched * TEMPERATURE_U,
    );
  }
}
