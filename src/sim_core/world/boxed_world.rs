use std::io;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::data::{ParticleConfig, SimulationConfig};

use crate::persistence::dto::world::WorldDTO;
use crate::persistence::dto::world::boxed::BoxedWorldDTO;

use crate::sim_core::world::WorldType;
use crate::sim_core::world::boxed_world::box_task::task_manager::TaskManager;
use crate::sim_core::world::boxed_world::computation_collector::ComputationCollector;
use crate::sim_core::world::boxed_world::history_manager::HistoryManager;
use crate::sim_core::world::boxed_world::integration_cache::IntegrationCache;
use crate::sim_core::world::boxed_world::integration_cache::integration_cache_builder::IntegrationCacheBuilder;
use crate::sim_core::world::boxed_world::persistance_reset::PersistanceReset;
use crate::sim_core::world::integration::{
  IntegrationAlgorithm, IntegrationAlgorithmState, new_integration_algorithm_state,
};
use crate::sim_core::world::saver::PartialWorldSaver;

pub mod box_container;
pub mod box_task;
mod computation_collector;
pub mod history_manager;
pub mod integration;
mod integration_cache;
pub mod persistance_reset;

pub struct BoxedWorld {
  config: SimulationConfig,
  persistance_reset: PersistanceReset,
  integration_cache_builder: IntegrationCacheBuilder,
  integration_cache: Option<Arc<IntegrationCache>>,
  computation_collector: Option<ComputationCollector>,
  task_manager: TaskManager,

  iteration: usize,

  // TODO: improve the integration_algorithm_logic
  integration_algorithm_state: IntegrationAlgorithmState,
}

impl BoxedWorld {
  pub fn with_config(config: SimulationConfig, particle_config: ParticleConfig) -> Self {
    let task_manager_config = match config.world_type {
      WorldType::BoxedWorld {
        task_manager_config,
      } => task_manager_config,
      _ => unreachable!(),
    };

    let num_of_atoms = particle_config.num_of_atoms;
    let history_manager = HistoryManager::with_config(config.clone(), particle_config);

    let box_container_config = *history_manager.box_container_config();

    let mut task_manager = TaskManager::new(
      task_manager_config.debug,
      config.clone(),
      box_container_config,
      task_manager_config.task_worker_multiplier,
    );

    task_manager.split_into_tasks_multiplier(&box_container_config);

    let initial_particles = history_manager
      .current_box_container()
      .all_particles_reset();
    let integration_cache_builder =
      IntegrationCacheBuilder::new(
        config.clone(),
        box_container_config, 
        initial_particles
      );

    let persistance_reset = PersistanceReset::new(
      history_manager,
      config.clone(),
      num_of_atoms,
      config.save_options.laamps_sampling.frame_iteration_count,
      config.save_options.energy_sampling.frame_iteration_count,
      PartialWorldSaver::new(config.save_options.clone()),
    );

    BoxedWorld {
      integration_algorithm_state: new_integration_algorithm_state(&config.integration_algorithm),
      config,
      persistance_reset,
      integration_cache_builder,
      integration_cache: None,
      computation_collector: None,
      task_manager,
      iteration: 0,
    }
  }

  pub fn save(&mut self) -> io::Result<()> {
    self
      .persistance_reset
      .save_full_snapshot_blocking(self.iteration)
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
        unimplemented!("semi implicit euler for boxed world")
      }
      IntegrationAlgorithm::VelocityVerlet => {
        unimplemented!("velocity verlet for boxed world!")
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

  pub fn get_size(&self) -> &Vector3<f64> {
    &self.config.world_size
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    let num_of_atoms = self
      .persistance_reset
      .history_manager()
      .get_particle_counts()
      .0;
    let lower_index: usize;
    if self.persistance_reset.number_of_resets() > 0 {
      lower_index = 1;
    } else {
      lower_index = 0;
    }

    WorldDTO::BoxedWorldDTO(BoxedWorldDTO {
      num_of_atoms,
      size: self.config.world_size,
      history: self
        .persistance_reset
        .history_manager()
        .to_transfer_struct(lower_index),
      integration_algorithm: self.config.integration_algorithm.clone(),

      num_of_world_iterations: self.iteration,
      number_of_resets: self.persistance_reset.number_of_resets(),
      max_iteration_till_reset: self.config.max_iteration_till_reset,

      laamps_frame_iteration_count: self.persistance_reset.laamps_frame_iteration_count(),
      energy_frame_iteration_count: self.persistance_reset.energy_frame_iteration_count(),
    })
  }

  pub fn get_particle_counts(&self) -> (usize, usize, usize) {
    self
      .persistance_reset
      .history_manager()
      .get_particle_counts()
  }
}
