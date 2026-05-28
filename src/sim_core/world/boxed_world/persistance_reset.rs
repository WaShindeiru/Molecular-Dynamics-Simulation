mod saver_handle;
mod saver_worker;

use std::fs;
use std::io;
use std::path::Path;
use std::sync::Arc;

use log::info;

use crate::data::types::AtomType;
use crate::data::{SimulationConfig, ValueUnits};
use crate::persistence::dto::world::boxed::BoxedWorldDTO;
use crate::persistence::dto::world::history::HistoryDTO;
use crate::persistence::json::particle_config::{ParticleConfigFile, ParticleInitialState};
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::history_manager::HistoryManager;
use crate::sim_core::world::saver::{PartialWorldSaver, SaveOptions};

use saver_handle::SaverHandle;

pub struct PersistanceReset {
  history_manager: HistoryManager,
  save_options: SaveOptions,
  config: SimulationConfig,
  num_of_atoms: usize,
  laamps_frame_iteration_count: usize,
  energy_frame_iteration_count: usize,
  reset_counter: usize,
  number_of_resets: usize,
  saver_handle: SaverHandle,
}

impl PersistanceReset {
  pub fn new(
    history_manager: HistoryManager,
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

  pub fn history_manager(&self) -> &HistoryManager {
    &self.history_manager
  }

  pub fn history_manager_mut(&mut self) -> &mut HistoryManager {
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

  fn build_boxed_snapshot_dto(&self, num_of_world_iterations: usize) -> BoxedWorldDTO {
    let lower_index = if self.number_of_resets > 0 { 1 } else { 0 };
    let (containers, thermostat_epsilon) = self.history_manager.clone_for_cache();
    let history = history_dto_from_snapshot(containers, thermostat_epsilon, lower_index);

    BoxedWorldDTO {
      num_of_atoms: self.num_of_atoms,
      size: self.config.world_size,
      history,
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
      if self.save_options.save {
        let dto = self.build_boxed_snapshot_dto(num_of_world_iterations);
        self.saver_handle.send_dto(dto)?;
      }

      info!(
        "Resetting boxed world, reset number: {}",
        self.number_of_resets
      );
      self.history_manager.reset_container();
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
      "Resetting boxed world, reset number: {}",
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
    let dto = self.build_boxed_snapshot_dto(num_of_world_iterations);
    self.saver_handle.send_dto_wait_result(dto)
  }

  pub fn save_final_particles(&self) -> io::Result<()> {
    if !self.save_options.save || !self.save_options.save_final_particles {
      return Ok(());
    }

    info!("Saving final particle states");

    let container = self.history_manager.current_box_container();

    let mut c_count = 0usize;
    let mut fe_count = 0usize;

    let particles: Vec<ParticleInitialState> = container
      .simulation_boxes()
      .iter()
      .flat_map(|sim_box| sim_box.particles().values())
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
    }
    .to_value_units(ValueUnits::Si);

    let dir = Path::new(&self.save_options.save_path).join("particles");
    fs::create_dir_all(&dir)?;

    let json = serde_json::to_string_pretty(&config_file)
      .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    fs::write(dir.join("final_particles.json"), json)
  }
}

impl Drop for PersistanceReset {
  fn drop(&mut self) {
    self.saver_handle.shutdown_inner();
  }
}

fn history_dto_from_snapshot(
  containers: Vec<Arc<BoxContainer<Arc<SimulationBox>>>>,
  thermostat_epsilon: Vec<f64>,
  lower_index: usize,
) -> HistoryDTO {
  let box_container = containers[lower_index..]
    .iter()
    .map(|c| c.to_transfer_struct())
    .collect();

  HistoryDTO {
    box_container,
    thermostat_epsilon,
  }
}
