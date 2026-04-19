use nalgebra::Vector3;
use crate::output::{EngineDTO};
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

  // pub fn new_from_atoms(atoms: Vec<Particle>, size: Vector3<f64>, time_step: f64,
  //                       num_of_iterations: usize,
  //                       max_iteration_till_reset: usize,
  //                       save_all_iterations: bool,
  //                       one_frame_duration: f64,
  //                       mut save_options: SaveOptions,
  //                       integration_algorithm: IntegrationAlgorithm,
  //                       world_type: WorldType,
  //                       edge_condition: EdgeCondition,
  // ) -> Self {
  // 
  //   Engine::new_from_atoms_with_path(atoms, size, time_step, num_of_iterations, max_iteration_till_reset, 
  //                                    save_all_iterations, one_frame_duration, save_options, 
  //                                    integration_algorithm, world_type, edge_condition
  //   )
  // }
  // 
  // pub fn new_from_atoms_with_path(atoms: Vec<Particle>, size: Vector3<f64>, time_step: f64,
  //                       num_of_iterations: usize,
  //                       max_iteration_till_reset: usize,
  //                       save_all_iterations: bool,
  //                       one_frame_duration: f64,
  //                       save_options: SaveOptions,
  //                       integration_algorithm: IntegrationAlgorithm,
  //                       world_type: WorldType, 
  //                       edge_condition: EdgeCondition
  // ) -> Self {
  //   let config = SimulationConfig::new(
  //     atoms,
  //     size,
  //     time_step,
  //     num_of_iterations,
  //     max_iteration_till_reset,
  //     save_all_iterations,
  //     one_frame_duration,
  //     save_options,
  //     integration_algorithm,
  //     world_type,
  //     edge_condition,
  //   );
  // 
  //   Self::from_config(config)
  // }

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
    let time_string = now.format("%Y-%m-%d_%H-%M-%S").to_string();

    let mut file = File::create(&format!("./{}/info.txt", self.save_options.save_path))?;
    writeln!(file, "Simulation date: {}", time_string)?;
    writeln!(file, "Number of iterations : {:?}", self.config.num_of_iterations)?;
    writeln!(file, "Time step: {:.3e} seconds", self.config.time_step * TIME_U)?;
    writeln!(file, "Simulation Time: {:?} seconds", self.simulation_time)?;
    writeln!(file, "integration type: {}", self.config.integration_algorithm)?;

    // Particle information
    let (total_particles, carbon_count, iron_count) = self.world.get_particle_counts();
    writeln!(file, "\n=== Particle Information ===")?;
    writeln!(file, "Total particles: {}", total_particles)?;
    writeln!(file, "Carbon (C) particles: {}", carbon_count)?;
    writeln!(file, "Iron (Fe) particles: {}", iron_count)?;

    // Gravity constant
    writeln!(file, "\n=== Physical Parameters ===")?;
    writeln!(file, "Gravity constant: {:.3e}", self.config.potential_gravity_max)?;

    // Integration-specific information
    writeln!(file, "\n=== Integration Algorithm Details ===")?;
    match &self.config.integration_algorithm {
      IntegrationAlgorithm::SemiImplicitEuler => {
        writeln!(file, "Algorithm: Semi-Implicit Euler")?;
      }
      IntegrationAlgorithm::VelocityVerlet => {
        writeln!(file, "Algorithm: Velocity Verlet")?;
      }
      IntegrationAlgorithm::NoseHooverVerlet { q_effective_mass, desired_temperature } => {
        writeln!(file, "Algorithm: Nose-Hoover Verlet")?;
        writeln!(file, "Q effective mass: {:.6e}", q_effective_mass)?;
        writeln!(file, "Temperature schedule:")?;
        for (idx, temp_info) in desired_temperature.iter().enumerate() {
          writeln!(file, "  Stage {}: {:.2} K, distance: {:?}",
            idx + 1,
            temp_info.desired_temperature,
            temp_info.distance
          )?;
        }
      }
    }

    writeln!(file, "\n=== World Configuration ===")?;
    writeln!(file, "World type: {}", self.world.get_world_info())?;
    match &self.world {
      World::SimpleWorld(simple_world) => {
        writeln!(file, "Simple world specific parameters:")?;
        writeln!(file, "  (No special parameters for SimpleWorld)")?;
      }
      World::BoxedWorld(boxed_world) => {
        writeln!(file, "Boxed world specific parameters:")?;
        writeln!(file, "  Box-based spatial partitioning active")?;
      }
    }

    Ok(())
  }
}
