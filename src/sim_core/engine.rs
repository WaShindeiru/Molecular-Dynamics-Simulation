use std::collections::HashMap;
use nalgebra::Vector3;
use crate::output::{AtomDTO, EngineDTO};
use crate::sim_core::world::{World, WorldType};

use std::{fs, io, time};
use std::error::Error;
use std::io::Write;
use std::time::{Duration, Instant};
use chrono::prelude::*;
use csv::Writer;
use log::info;
use crate::particle::Particle;
use crate::sim_core::world::integration::{IntegrationAlgorithm, IntegrationAlgorithmParams, validate_integration_params};
use crate::sim_core::world::saver::SaveOptions;

pub struct Engine {
  world: World,
  time_step: f64,
  current_time: f64,
  current_iteration: usize,
  num_of_iterations: usize,
  simulation_time: Duration,
  integration_algorithm: IntegrationAlgorithm,
  
  save_all_iterations: bool,
  one_frame_duration: f64,
  frame_iteration_count: usize,
  max_iteration_till_reset: usize,

  save_options: SaveOptions,
}

impl Engine {
  pub fn new_from_atoms(atoms: Vec<Particle>, size: Vector3<f64>, time_step: f64,
                        num_of_iterations: usize,
                        max_iteration_till_reset: usize,
                        save_all_iterations: bool,
                        one_frame_duration: f64,
                        mut save_options: SaveOptions,
                        integration_algorithm: IntegrationAlgorithm,
                        world_type: WorldType,
  ) -> Self {
    
    let now: DateTime<Local> = Local::now();
    let time_string = now.format("%Y-%m-%d_%H-%M-%S").to_string();
    let save_path = "../output/".to_string() + &*time_string;
    save_options.save_path = save_path;

    Engine::new_from_atoms_with_path(atoms, size, time_step, num_of_iterations, max_iteration_till_reset, 
                                     save_all_iterations, one_frame_duration, save_options, 
                                     integration_algorithm, world_type
    )
  }

  pub fn new_from_atoms_with_path(atoms: Vec<Particle>, size: Vector3<f64>, time_step: f64,
                        num_of_iterations: usize,
                        max_iteration_till_reset: usize,
                        save_all_iterations: bool,
                        one_frame_duration: f64,
                        save_options: SaveOptions,
                        integration_algorithm: IntegrationAlgorithm,
                        world_type: WorldType,
  ) -> Self {
    let mut frame_iteration_count = 1;

    if !save_all_iterations {
      frame_iteration_count = (one_frame_duration / time_step) as usize;
    }

    let world = World::new_from_atoms(atoms, size, max_iteration_till_reset, frame_iteration_count,
                                      integration_algorithm.clone(), save_options.clone(), world_type);

    Engine {
      world,
      time_step,
      current_time: 0.0,
      current_iteration: 0,
      num_of_iterations,
      simulation_time: Duration::ZERO,
      integration_algorithm,

      save_all_iterations,
      one_frame_duration,
      frame_iteration_count,
      max_iteration_till_reset,

      save_options,
    }
  }

  pub fn run(&mut self, params: &IntegrationAlgorithmParams, time_step: f64) {
    assert!(
      validate_integration_params(&self.integration_algorithm, params),
      "IntegrationAlgorithmParams {:?} does not match IntegrationAlgorithm {}",
      params, self.integration_algorithm
    );

    let start = Instant::now();
    let spinner = ['|', '/', '-', '\\'];
    let mut counter = 0;
    
    for i in 0..self.num_of_iterations {
      if i % 100 == 0 {
        let frame = counter % spinner.len();
        print!("\rProgress: {}/{} {}", i, self.num_of_iterations, spinner[frame]);
        io::stdout().flush().unwrap();
        counter += 1;
      }

      self.world.update(params, time_step, self.current_iteration + 1);

      self.current_iteration += 1;
      self.current_time += self.time_step;
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
      time_step: self.time_step,
      num_of_iterations: self.num_of_iterations,
    }
  }

  pub fn save(&mut self) -> io::Result<()> {
    self.world.save()?;

    let now: DateTime<Local> = Local::now();
    let time_string = now.format("%Y-%m-%d_%H-%M-%S").to_string();

    let mut wtr = Writer::from_path(&format!("./{}/info.txt", self.save_options.save_path))?;
    wtr.write_record(&["Simulation date: ", &format!("{}", time_string)])?;
    wtr.write_record(&["Number of iterations : ", &format!("{:.2?}", self.num_of_iterations)])?;
    wtr.write_record(&["Time step: ", &format!("{:.2?} seconds", self.time_step)])?;
    wtr.write_record(&["Simulation Time: ", &format!("{:.2?} seconds", self.simulation_time)])?;
    wtr.write_record(&["integration type: ", &format!("{}", self.integration_algorithm)])?;

    match self.world {
      World::SimpleWorld(_) => wtr.write_record(&["world type: Simple World", ""])?,
      World::BoxedWorld(_) => wtr.write_record(&["world type: Boxed World", ""])?,
    };

    wtr.flush()?;

    Ok(())
  }
}
