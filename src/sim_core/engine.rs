use std::collections::HashMap;
use nalgebra::Vector3;
use crate::output::{AtomDTO, EngineDTO};
use crate::sim_core::world::World;

use std::{fs, io, time};
use std::io::Write;
use std::time::{Duration, Instant};
use chrono::prelude::*;
use csv::Writer;
use log::info;
use crate::particle::Particle;
use crate::sim_core::world::integration::IntegrationAlgorithm;

pub struct Engine {
  world: World,
  time_step: f64,
  current_time: f64,
  current_iteration: usize,
  num_of_iterations: usize,
  simulation_time: Duration,
  
  save_all_iterations: bool,
  one_frame_duration: f64,
  frame_iteration_count: usize,
  max_iteration_till_reset: usize,

  save: bool,
  save_path: String,
  save_laamps: bool,
  save_verbose: bool,
}

impl Engine {
  // pub fn new(world: World, time_step: f64, num_of_iterations: usize) -> Self {
  //   Engine {
  //     world,
  //     time_step,
  //     current_time: 0.0,
  //     current_iteration: 0,
  //     num_of_iterations,
  //     simulation_time: Duration::ZERO
  //   }
  // }

  pub fn new_from_atoms(atoms: Vec<Particle>, size: Vector3<f64>, time_step: f64,
                        num_of_iterations: usize,
                        max_iteration_till_reset: usize,
                        save: bool,
                        save_laamps: bool,
                        save_verbose: bool,
                        save_all_iterations: bool,
                        one_frame_duration: f64, ) -> Self {
    let mut frame_iteration_count = 1;
    
    if !save_all_iterations {
      frame_iteration_count = (one_frame_duration / time_step) as usize;
    }

    let now: DateTime<Local> = Local::now();
    let time_string = now.format("%Y-%m-%d_%H-%M-%S").to_string();
    let save_path = "../output/".to_string() + &*time_string;

    let world = World::new_from_atoms(atoms, size, max_iteration_till_reset, frame_iteration_count,
                                            save, save_path.clone(), save_laamps, save_verbose);

    Engine {
      world,
      time_step,
      current_time: 0.0,
      current_iteration: 0,
      num_of_iterations,
      simulation_time: Duration::ZERO,
      save_all_iterations,
      one_frame_duration,
      frame_iteration_count,
      max_iteration_till_reset,
      save,
      save_path,
      save_laamps,
      save_verbose,
    }
  }
  
  pub fn run(&mut self, params: &IntegrationAlgorithm, time_step: f64) {
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

    if self.save {
      let use_thermostat: bool;
      match params {
        IntegrationAlgorithm::NoseHooverVerlet {..} => {
          use_thermostat = true;
        }
        _ => {
          use_thermostat = false;
        }
      }

      self.save(use_thermostat).unwrap()
    }
  }

  pub fn to_transfer_struct(&self) -> EngineDTO {
    EngineDTO {
      world: self.world.to_transfer_struct(),
      time_step: self.time_step,
      num_of_iterations: self.num_of_iterations,
    }
  }

  pub fn save(&mut self, use_thermostat: bool) -> io::Result<()> {
    self.world.save(use_thermostat)?;

    let now: DateTime<Local> = Local::now();
    let time_string = now.format("%Y-%m-%d_%H-%M-%S").to_string();

    let mut wtr = Writer::from_path(&format!("./{}/info.txt", self.save_path))?;
    wtr.write_record(&["Simulation date: ", &format!("{}", time_string)])?;
    wtr.write_record(&["Number of iterations : ", &format!("{:.2?}", self.num_of_iterations)])?;
    wtr.write_record(&["Time step: ", &format!("{:.2?} seconds", self.time_step)])?;
    wtr.write_record(&["Simulation Time: ", &format!("{:.2?} seconds", self.simulation_time)])?;
    wtr.write_record(&["use thermostat: ", &format!("{}", use_thermostat)])?;
    wtr.flush()?;

    Ok(())
  }
}
