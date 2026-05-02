use nalgebra::Vector3;
use crate::output::engine::EngineDTO;
use crate::sim_core::world::{World, WorldType};
use crate::data::SimulationConfig;

use std::io;
use std::io::Write;
use std::time::{Duration, Instant};
use chrono::prelude::*;
use std::fs::File;
use log::info;
use crate::data::units::TIME_U;
use crate::particle::Particle;
use crate::sim_core::world::integration::{IntegrationAlgorithm};
use crate::sim_core::world::saver::SaveOptions;
use crate::sim_core::world::boundary_constraint::EdgeCondition;

pub struct Engine {
  config: SimulationConfig,
  world: World,

  current_time: f64,
  current_iteration: usize,
  simulation_time: Duration,

  save_options: SaveOptions,
}

impl Engine {
  pub fn from_config(config: SimulationConfig) -> Self {
    let world = World::from_config(config.clone());

    Engine {
      config: config.clone(),
      world,
      current_time: 0.0,
      current_iteration: 0,
      simulation_time: Duration::ZERO,

      save_options: config.save_options.clone(),
    }
  }

  pub fn run(&mut self, time_step: f64) {
    info!("Starting simulation...");
    info!("Will save output to {}", self.save_options.save_path);

    let start = Instant::now();
    let spinner = ['|', '/', '-', '\\'];
    let mut counter = 0;
    
    for i in 0..self.config.num_of_iterations {
      if i % 100 == 0 {
        let frame = counter % spinner.len();
        print!("\rProgress: {}/{} {}", i, self.config.num_of_iterations, spinner[frame]);
        io::stdout().flush().unwrap();
        counter += 1;
      }

      self.world.update(&self.config.integration_algorithm, time_step, self.current_iteration + 1);

      self.current_iteration += 1;
      self.current_time += self.config.time_step;
    }

    print!("\n");
    io::stdout().flush().unwrap();

    self.simulation_time = start.elapsed();

    info!("Simulation completed in {:.2?} seconds.", self.simulation_time);

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

    let now: DateTime<Local> = Local::now();
    let time_string: String = now.format("%Y-%m-%d_%H-%M-%S").to_string();

    let mut file = File::create(&format!("./{}/info.txt", self.save_options.save_path))?;
    writeln!(file, "Simulation date: {}", time_string)?;
    writeln!(file, "Number of iterations : {:?}", self.config.num_of_iterations)?;
    writeln!(file, "Simulation Time: {:?} seconds", self.simulation_time)?;

    self.config.to_json_file(&format!("./{}/parameters.json", self.save_options.save_path))?;

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
