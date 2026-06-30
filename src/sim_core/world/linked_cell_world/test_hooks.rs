use crate::sim_core::world::linked_cell_world::LinkedCellContainer;
use crate::sim_core::world::linked_cell_world::LinkedCellWorld;
use crate::sim_core::world::linked_cell_world::integration::verlet_nose_hoover::thermostat::compute_new_thermostat_epsilon;
use crate::sim_core::world::thermostat::{
    IntegrationAlgorithm, IntegrationAlgorithmState, IntegrationStateUpdateResponse,
};

impl LinkedCellWorld {
    /// Instrumented variant of `update_verlet_nose_hoover` for testing.
    /// Calls the three callbacks at checkpoints without altering the update logic:
    ///   `on_before`             — called with the container at the start of the step
    ///   `on_after_half_velocity` — called with the container after the half-velocity step
    ///   `on_after_update`       — called with the container at the end of the step
    pub fn update_verlet_nose_hoover_instrumented<F1, F2, F3>(
        &mut self,
        next_iteration: usize,
        on_before: F1,
        on_after_half_velocity: F2,
        on_after_update: F3,
    ) where
        F1: FnOnce(&LinkedCellContainer),
        F2: FnOnce(&LinkedCellContainer),
        F3: FnOnce(&LinkedCellContainer),
    {
        let current_thermostat_epsilon = self
            .persistance_reset
            .history_manager()
            .current_thermostat_epsilon();
        let current_container = self
            .persistance_reset
            .history_manager()
            .current_container();

        on_before(&current_container);

        let temperature_index =
            if let IntegrationAlgorithmState::NoseHooverVerlet { temperature_index, .. } =
                self.integration_algorithm_state
            {
                temperature_index
            } else {
                panic!("Expected NoseHooverVerlet integration state!")
            };

        let (current_desired_temperature, q_effective_mass) =
            if let IntegrationAlgorithm::NoseHooverVerlet {
                desired_temperature,
                q_effective_mass: q,
            } = &self.config.integration_algorithm
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

        on_after_half_velocity(&integration_cache.local_container);

        let new_thermostat_epsilon = compute_new_thermostat_epsilon(
            current_thermostat_epsilon,
            &integration_cache.half_velocity_cache,
            integration_cache
                .local_container
                .particles()
                .iter()
                .filter_map(|p| p.as_ref()),
            self.config.time_step,
            q_effective_mass,
            current_desired_temperature,
        );

        self.persistance_reset
            .history_manager_mut()
            .add_thermostat_epsilon(new_thermostat_epsilon);

        let mut computation_collector = self.task_manager.force_step(integration_cache);

        computation_collector.apply_gravity(next_iteration);
        computation_collector.set_velocity(new_thermostat_epsilon);

        let current_custom_velocities = self
            .velocity_manager
            .compute_velocities_for_iteration(next_iteration);
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
            IntegrationStateUpdateResponse::NoseHooverVerlet { .. } => {}
            _ => panic!("Wrong result type"),
        }

        self.persistance_reset
            .history_manager_mut()
            .push_container(computation_collector.build());

        self.iteration += 1;
        assert_eq!(self.iteration, next_iteration);

        let final_container = self
            .persistance_reset
            .history_manager()
            .current_container();

        on_after_update(&final_container);
    }
}
