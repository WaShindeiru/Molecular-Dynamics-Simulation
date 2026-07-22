pub mod computation_collector;
pub mod history_manager;
pub mod integration;
pub mod integration_cache;
pub mod linked_cell_container;
pub mod linked_cell_task;
pub mod persistance_reset;

pub mod test_hooks;

pub use history_manager::LinkedCellHistoryManager;
pub use linked_cell_container::LinkedCellContainerOld;
pub use persistance_reset::LinkedCellPersistanceReset;

use std::io;
use std::path::Path;

use crate::data::{ParticleConfig, SimulationConfig};
use crate::particle::Particle;
use crate::persistence::dto::world::WorldDTO;
use crate::persistence::dto::world::boxed::BoxedWorldDTO;
use crate::persistence::json::control_velocity_manager_file::ControlVelocityManagerFile;
use crate::persistence::json::velocity_manager_file::VelocityManagerFile;
use crate::sim_core::world::WorldType;
use crate::sim_core::world::boxed_world::velocity_manager::VelocityManager;
use crate::sim_core::world::linked_cell_world::linked_cell_task::task_manager::TaskManager;
use crate::sim_core::world::saver::PartialWorldSaver;
use crate::sim_core::world::thermostat::{IntegrationAlgorithm, IntegrationAlgorithmState, new_integration_algorithm_state};

pub struct LinkedCellWorld {
  pub config: SimulationConfig,
  pub persistance_reset: LinkedCellPersistanceReset,
  pub task_manager: TaskManager,
  pub velocity_manager: VelocityManager,
  pub iteration: usize,
  pub integration_algorithm_state: IntegrationAlgorithmState,
}

impl LinkedCellWorld {
  pub fn with_config(config: SimulationConfig, mut particle_config: ParticleConfig) -> Self {
    let task_manager_config = match config.world_type {
      WorldType::LinkedCellWorld { task_manager_config } => task_manager_config,
      _ => unreachable!(),
    };

    assert!(
      !particle_config
        .atoms
        .iter()
        .any(|atom| matches!(atom, Particle::VelocityControlledParticle(_))),
      "VelocityControlledParticle is only supported in OptimizedWorld"
    );

    let mut velocity_manager = VelocityManager::from_config(&particle_config);
    let initial_velocities = velocity_manager.compute_velocities_for_iteration(0);
    for &(particle_id, vel) in initial_velocities {
      if let Particle::CustomVelocityAtom(p) = &mut particle_config.atoms[particle_id] {
        p.set_velocity(vel);
      }
    }

    let velocity_managers_file: Vec<VelocityManagerFile> = particle_config
      .velocity_schedules
      .iter()
      .map(VelocityManagerFile::from_schedule)
      .collect();
    let control_velocity_managers_file: Vec<ControlVelocityManagerFile> = particle_config
      .control_velocity_schedules
      .iter()
      .map(ControlVelocityManagerFile::from_schedule)
      .collect();

    let num_of_atoms = particle_config.num_of_atoms;
    let history_manager = LinkedCellHistoryManager::with_config(config.clone(), particle_config);

    let current_container = history_manager.current_container();

    let mut task_manager = TaskManager::new(
      task_manager_config.debug,
      config.clone(),
      task_manager_config.task_worker_multiplier,
      task_manager_config.split,
      num_of_atoms,
    );

    task_manager.split_into_tasks_multiplier(&current_container);

    let persistance_reset = LinkedCellPersistanceReset::new(
      history_manager,
      config.clone(),
      num_of_atoms,
      config.save_options.laamps_sampling.frame_iteration_count,
      config.save_options.energy_sampling.frame_iteration_count,
      PartialWorldSaver::new(
        config.save_options.clone(),
        velocity_managers_file.clone(),
        control_velocity_managers_file.clone(),
      ),
      velocity_managers_file,
      control_velocity_managers_file,
    );

    LinkedCellWorld {
      integration_algorithm_state: new_integration_algorithm_state(&config.integration_algorithm),
      config,
      persistance_reset,
      task_manager,
      velocity_manager,
      iteration: 0,
    }
  }

  pub fn save(&mut self) -> io::Result<()> {
    if self.config.save_options.save_final_particles {
      self.persistance_reset.save_final_particles()?;
    }

    self.integration_algorithm_state.save_temperature_particles(
      &self.config.integration_algorithm,
      &Path::new(&self.config.save_options.save_path),
      self.persistance_reset.velocity_managers_file(),
      self.persistance_reset.control_velocity_managers_file(),
    )?;

    self.persistance_reset.save_full_snapshot_blocking(self.iteration)
  }

  pub fn update(
    &mut self,
    algorithm: &IntegrationAlgorithm,
    _time_step: f64,
    next_iteration: usize,
  ) -> io::Result<()> {
    self.persistance_reset.begin_update_step(self.iteration)?;

    match algorithm {
      IntegrationAlgorithm::SemiImplicitEuler => {
        unimplemented!("semi implicit euler for linked cell world")
      }
      IntegrationAlgorithm::VelocityVerlet => {
        self.update_velocity_verlet(next_iteration);
      }
      IntegrationAlgorithm::NoseHooverVerlet { .. } => {
        self.update_verlet_nose_hoover(next_iteration);
      }
    }

    self.persistance_reset.end_update_step();

    Ok(())
  }

  pub fn reset_world(&mut self) {
    self.persistance_reset.reset_world_without_save();
  }

  pub fn get_size(&self) -> &nalgebra::Vector3<f64> {
    &self.config.world_size
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    let num_of_atoms = self.persistance_reset.history_manager().get_particle_counts().0;
    let lower_index = if self.persistance_reset.number_of_resets() > 0 { 1 } else { 0 };

    WorldDTO::BoxedWorldDTO(BoxedWorldDTO {
      num_of_atoms,
      size: self.config.world_size,
      history: self.persistance_reset.history_manager().to_transfer_struct(lower_index),
      integration_algorithm: self.config.integration_algorithm.clone(),
      num_of_world_iterations: self.iteration,
      number_of_resets: self.persistance_reset.number_of_resets(),
      max_iteration_till_reset: self.config.max_iteration_till_reset,
      laamps_frame_iteration_count: self.persistance_reset.laamps_frame_iteration_count(),
      energy_frame_iteration_count: self.persistance_reset.energy_frame_iteration_count(),
    })
  }

  pub fn get_particle_counts(&self) -> (usize, usize, usize) {
    self.persistance_reset.history_manager().get_particle_counts()
  }
}
