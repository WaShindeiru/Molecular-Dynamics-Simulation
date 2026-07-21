mod saver_handle;
mod saver_worker;

use std::fs;
use std::io;
use std::path::Path;

use log::info;

use crate::data::types::AtomType;
use crate::data::{SimulationConfig, ValueUnits};
use crate::perf_log;
use crate::persistence::dto::world::boxed::BoxedWorldDTOWithoutHistory;
use crate::persistence::json::particle_config::{ParticleConfigFile, ParticleInitialState};
use crate::sim_core::world::optimized_world::history_manager::OptimizedHistoryManager;
use crate::sim_core::world::saver::{PartialWorldSaver, SaveOptions};

use saver_handle::SaverHandle;

pub struct OptimizedPersistanceReset {
  history_manager: OptimizedHistoryManager,
  save_options: SaveOptions,
  config: SimulationConfig,
  num_of_atoms: usize,
  laamps_frame_iteration_count: usize,
  energy_frame_iteration_count: usize,
  reset_counter: usize,
  number_of_resets: usize,
  saver_handle: SaverHandle,
}

impl OptimizedPersistanceReset {
  pub fn new(
    history_manager: OptimizedHistoryManager,
    config: SimulationConfig,
    num_of_atoms: usize,
    laamps_frame_iteration_count: usize,
    energy_frame_iteration_count: usize,
    saver: PartialWorldSaver,
  ) -> Self {
    let save_enabled = config.save_options.save;
    let save_options = config.save_options.clone();
    Self {
      history_manager,
      save_options,
      config,
      num_of_atoms,
      laamps_frame_iteration_count,
      energy_frame_iteration_count,
      reset_counter: 1,
      number_of_resets: 0,
      saver_handle: SaverHandle::spawn(saver, save_enabled),
    }
  }

  pub fn history_manager(&self) -> &OptimizedHistoryManager {
    &self.history_manager
  }

  pub fn history_manager_mut(&mut self) -> &mut OptimizedHistoryManager {
    &mut self.history_manager
  }

  pub fn number_of_resets(&self) -> usize {
    self.number_of_resets
  }

  pub fn laamps_frame_iteration_count(&self) -> usize {
    self.laamps_frame_iteration_count
  }

  pub fn energy_frame_iteration_count(&self) -> usize {
    self.energy_frame_iteration_count
  }

  fn build_partial_dto(&self, num_of_world_iterations: usize) -> BoxedWorldDTOWithoutHistory {
    BoxedWorldDTOWithoutHistory {
      num_of_atoms: self.num_of_atoms,
      size: self.config.world_size,
      integration_algorithm: self.config.integration_algorithm.clone(),
      num_of_world_iterations,
      number_of_resets: self.number_of_resets,
      max_iteration_till_reset: self.config.max_iteration_till_reset,
      laamps_frame_iteration_count: self.laamps_frame_iteration_count,
      energy_frame_iteration_count: self.energy_frame_iteration_count,
    }
  }

  pub fn begin_update_step(&mut self, num_of_world_iterations: usize) -> io::Result<()> {
    assert!(self.reset_counter <= self.config.max_iteration_till_reset);

    self.saver_handle.drain_finished_results()?;

    if self.reset_counter == self.config.max_iteration_till_reset {
      perf_log!("Performing container reset");
      if self.save_options.save {
        let partial = self.build_partial_dto(num_of_world_iterations);
        let fresh = self.history_manager.reset_clone();
        let old = std::mem::replace(&mut self.history_manager, fresh);
        self.saver_handle.send_msg((partial, old))?;
      } else {
        self.history_manager.reset_container();
      }

      info!(
        "Resetting optimized world, reset number: {}",
        self.number_of_resets
      );
      self.number_of_resets += 1;
      self.reset_counter = 0;
    }

    Ok(())
  }

  pub fn end_update_step(&mut self) {
    self.reset_counter += 1;
  }

  pub fn reset_world_without_save(&mut self) {
    info!(
      "Resetting optimized world, reset number: {}",
      self.number_of_resets
    );
    self.history_manager.reset_container();
    self.number_of_resets += 1;
    self.reset_counter = 0;
  }

  pub fn save_full_snapshot_blocking(&mut self, num_of_world_iterations: usize) -> io::Result<()> {
    if !self.save_options.save {
      return Ok(());
    }
    self.saver_handle.drain_finished_results()?;
    let partial = self.build_partial_dto(num_of_world_iterations);
    let fresh = self.history_manager.reset_clone();
    let old = std::mem::replace(&mut self.history_manager, fresh);
    self.saver_handle.send_msg_wait_result((partial, old))
  }

  pub fn save_final_particles(&self) -> io::Result<()> {
    if !self.save_options.save || !self.save_options.save_final_particles {
      return Ok(());
    }

    info!("Saving final particle states");

    let container = self.history_manager.current_container();

    let mut c_count = 0usize;
    let mut fe_count = 0usize;

    let particles: Vec<ParticleInitialState> = container
      .particles()
      .iter()
      .map(|particle| {
        match particle.get_type() {
          AtomType::C | AtomType::C_nanotube => c_count += 1,
          AtomType::Fe => fe_count += 1,
        }
        ParticleInitialState::from_runtime(particle)
      })
      .collect();

    let config_file = ParticleConfigFile {
      value_units: ValueUnits::Unitless,
      particles,
      num_of_atoms: c_count + fe_count,
      num_of_carbon_atoms: c_count,
      num_of_iron_atoms: fe_count,
      velocity_managers: vec![],
      control_velocity_managers: vec![],
    }
    .to_value_units(ValueUnits::Si);

    let dir = Path::new(&self.save_options.save_path).join("particles");
    fs::create_dir_all(&dir)?;

    let json = serde_json::to_string_pretty(&config_file)
      .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    fs::write(dir.join("final_particles.json"), json)
  }
}

impl Drop for OptimizedPersistanceReset {
  fn drop(&mut self) {
    self.saver_handle.shutdown_inner();
  }
}
