use crate::data::{ConfigAll, ParticleConfig, SimulationConfig, ValueUnits};
use crate::persistence::dto::engine::EngineDTO;
use crate::persistence::json::SimulationConfigFile;
use crate::persistence::json::particle_config::particle_config_to_initial_json_file;
use crate::sim_core::world::World;

use crate::sim_core::world::saver::SaveOptions;
use chrono::prelude::*;
use log::info;
use std::fs::File;
use std::io;
use std::io::Write;
use std::path::Path;
use std::time::{Duration, Instant};

pub struct Engine {
  config: SimulationConfig,
  particle_config: ParticleConfig,
  world: World,

  current_time: f64,
  current_iteration: usize,
  simulation_time: Duration,

  save_options: SaveOptions,
}

impl Engine {
  pub fn from_configs(config: SimulationConfig, particle_config: ParticleConfig) -> Self {
    let world = World::from_configs(config.clone(), particle_config.clone());

    Engine {
      config: config.clone(),
      particle_config,
      world,
      current_time: 0.0,
      current_iteration: 0,
      simulation_time: Duration::ZERO,

      save_options: config.save_options.clone(),
    }
  }

  pub fn from_config_all(config_all: ConfigAll) -> Self {
    Self::from_configs(config_all.simulation_config, config_all.particle_config)
  }

  pub fn run(&mut self) {
    info!("Starting simulation...");
    info!("Will save output to {}", self.save_options.save_path);

    let start = Instant::now();
    let spinner = ['|', '/', '-', '\\'];
    let mut counter = 0;

    for i in 0..self.config.num_of_iterations {
      if i % 100 == 0 {
        let frame = counter % spinner.len();
        print!(
          "\rProgress: {}/{} {}",
          i, self.config.num_of_iterations, spinner[frame]
        );
        io::stdout().flush().unwrap();
        counter += 1;
      }

      self
        .world
        .update(
          &self.config.integration_algorithm,
          self.config.time_step,
          self.current_iteration + 1,
        )
        .unwrap();

      self.current_iteration += 1;
      self.current_time += self.config.time_step;
    }

    print!("\n");
    io::stdout().flush().unwrap();

    self.simulation_time = start.elapsed();

    info!(
      "Simulation completed in {:.2?} seconds.",
      self.simulation_time
    );

    if self.save_options.save {
      self.save().unwrap()
    }
  }

  pub fn to_transfer_struct(&self) -> EngineDTO {
    EngineDTO {
      world: self.world.to_transfer_struct(),
      time_step: self.config.time_step,
      num_of_iterations: self.config.num_of_iterations,
    }
  }

  pub fn save(&mut self) -> io::Result<()> {
    self.world.save()?;
    let save_dir = Path::new(&self.save_options.save_path);

    let now: DateTime<Local> = Local::now();
    let time_string: String = now.format("%Y-%m-%d_%H-%M-%S").to_string();

    let mut file = File::create(save_dir.join("info.txt"))?;
    writeln!(file, "Simulation date: {}", time_string)?;
    writeln!(
      file,
      "Number of iterations : {:?}",
      self.config.num_of_iterations
    )?;
    writeln!(file, "Simulation Time: {:?} seconds", self.simulation_time)?;

    let parameters_path = save_dir.join("parameters.json");
    SimulationConfigFile::from_runtime(&self.config, ValueUnits::Si)
      .to_json_file(parameters_path.to_string_lossy().as_ref())?;
    let particles_initial_path = save_dir.join("particles_initial.json");
    particle_config_to_initial_json_file(
      &self.particle_config,
      particles_initial_path.to_string_lossy().as_ref(),
    )?;

    // writeln!(file, "\n=== World Configuration ===")?;
    // writeln!(file, "World type: {}", self.world.get_world_info())?;
    // match &self.world {
    //   World::SimpleWorld(simple_world) => {
    //     writeln!(file, "Simple world specific parameters:")?;
    //     writeln!(file, "  (No special parameters for SimpleWorld)")?;
    //   }
    //   World::BoxedWorld(boxed_world) => {
    //     writeln!(file, "Boxed world specific parameters:")?;
    //     writeln!(file, "  Box-based spatial partitioning active")?;
    //   }
    // }

    Ok(())
  }
}
