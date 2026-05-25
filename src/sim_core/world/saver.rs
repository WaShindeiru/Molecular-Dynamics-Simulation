use std::io;
use std::path::Path;

use crate::persistence::dto::atom::AtomDTO;
use crate::persistence::dto::world::WorldDTO;
use crate::persistence::json::particle_config::particle_type_file::ParticleTypeFile;
use crate::persistence::json::velocity_particle::VelocityParticleFile;
use crate::sim_core::world::velocity_heap::{VelocityHeap, VelocityParticle};
use std::fs;

mod boxed_world;
mod simple_world;

#[derive(Clone, Copy)]
pub struct FrameSamplingConfig {
  pub one_frame_duration: f64,
  pub frame_iteration_count: usize,
}

impl FrameSamplingConfig {
  /// Nearest integer to `one_frame_duration / time_step` (half-way cases round away from zero).
  /// Always at least `1` for finite positive ratios.
  pub fn iterations_from_duration(one_frame_duration: f64, time_step: f64) -> usize {
    let ratio = one_frame_duration / time_step;
    if !ratio.is_finite() || ratio <= 0.0 {
      return 1;
    }
    let rounded = ratio.round();
    if !rounded.is_finite() || rounded <= 0.0 {
      return 1;
    }
    let capped = rounded.min(usize::MAX as f64);
    (capped as usize).max(1)
  }

  pub fn from_duration(time_step: f64, one_frame_duration: f64, save_all_iterations: bool) -> Self {
    let frame_iteration_count = if save_all_iterations {
      1
    } else {
      Self::iterations_from_duration(one_frame_duration, time_step)
    };

    FrameSamplingConfig {
      one_frame_duration,
      frame_iteration_count,
    }
  }
}

impl Default for FrameSamplingConfig {
  fn default() -> Self {
    FrameSamplingConfig {
      one_frame_duration: 1e-16,
      frame_iteration_count: 1,
    }
  }
}

#[derive(Clone)]
pub struct SaveOptions {
  pub save: bool,
  pub save_path: String,
  pub save_laamps: bool,
  pub save_verbose: bool,
  pub save_all_iterations_laamps: bool,
  pub save_all_iterations_energy: bool,
  pub laamps_sampling: FrameSamplingConfig,
  pub energy_sampling: FrameSamplingConfig,
  pub velocity_particles_num: usize,
  pub save_final_particles: bool,
}

impl Default for SaveOptions {
  fn default() -> Self {
    SaveOptions {
      save: false,
      save_path: String::from("output"),
      save_laamps: false,
      save_verbose: false,
      save_all_iterations_laamps: false,
      save_all_iterations_energy: false,
      laamps_sampling: FrameSamplingConfig::default(),
      energy_sampling: FrameSamplingConfig::default(),
      velocity_particles_num: 100,
      save_final_particles: false,
    }
  }
}

pub struct PartialWorldSaver {
  save_options: SaveOptions,

  thermostat_work_total: f64,
  laamps_frame_iteration_count_current_iteration: usize,
  energy_frame_iteration_count_current_iteration: usize,

  velocity_heap: VelocityHeap,
}

impl PartialWorldSaver {
  pub fn new(save_options: SaveOptions) -> Self {
    let velocity_heap = VelocityHeap::new(save_options.velocity_particles_num);
    PartialWorldSaver {
      thermostat_work_total: 0.,
      laamps_frame_iteration_count_current_iteration: 0,
      energy_frame_iteration_count_current_iteration: 0,
      velocity_heap,
      save_options,
    }
  }

  pub fn handle_velocity_particles(&mut self, atoms: impl Iterator<Item = AtomDTO>) {
    for atom in atoms {
      self.velocity_heap.try_insert(atom);
    }
  }

  pub fn persist_velocity_particles(&mut self) -> io::Result<()> {
    let dir = Path::new(&self.save_options.save_path).join("particles");
    fs::create_dir_all(&dir)?;

    let particles: Vec<VelocityParticleFile> = self
      .velocity_heap
      .drain_descending()
      .map(|vp| velocity_particle_to_file(vp).to_value_units(crate::data::ValueUnits::Si))
      .collect();

    let json = serde_json::to_string_pretty(&particles)
      .map_err(|e| io::Error::new(io::ErrorKind::InvalidData, e))?;

    fs::write(dir.join("velocity_particles.json"), json)
  }

  pub fn persist(&mut self, world: &WorldDTO) -> io::Result<()> {
    if self.save_options.save {
      match world {
        WorldDTO::SimpleWorldDTO(simple_world) => self.persist_simple_world(simple_world),
        WorldDTO::BoxedWorldDTO(boxed_world) => self.persist_boxed_world(boxed_world),
      }
    } else {
      Ok(())
    }
  }
}

fn velocity_particle_to_file(vp: VelocityParticle) -> VelocityParticleFile {
  VelocityParticleFile::new(
    vp.id,
    vp.iteration,
    vp.type_,
    ParticleTypeFile::from(vp.kind),
    vp.position,
    vp.velocity,
  )
}
