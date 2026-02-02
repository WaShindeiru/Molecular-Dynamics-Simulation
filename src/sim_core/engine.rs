use std::collections::HashMap;
use nalgebra::Vector3;
use crate::output::EngineDTO;
use crate::sim_core::world::World;

use std::{fs, io, time};
use std::io::Write;
use std::time::{Duration, Instant};
use chrono::prelude::*;
use csv::Writer;
use log::info;
use crate::particle::Particle;

pub struct Engine {
  world: World,
  time_step: f64,
  current_time: f64,
  current_iteration: usize,
  num_of_iterations: usize,
  simulation_time: time::Duration,
}

impl Engine {
  pub fn new(world: World, time_step: f64, num_of_iterations: usize) -> Self {
    Engine {
      world,
      time_step,
      current_time: 0.0,
      current_iteration: 0,
      num_of_iterations,
      simulation_time: Duration::ZERO
    }
  }

  pub fn new_from_atoms(atoms: Vec<Particle>, size: Vector3<f64>, time_step: f64, num_of_iterations: usize) -> Self {
    let world = World::new_from_atoms(atoms, size);

    Engine {
      world,
      time_step,
      current_time: 0.0,
      current_iteration: 0,
      num_of_iterations,
      simulation_time: Duration::ZERO
    }
  }
  
  pub fn run(&mut self, save: bool) {
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
      self.world.update_verlet(self.time_step, self.current_iteration + 1);
      self.current_iteration += 1;
      self.current_time += self.time_step;
    }

    print!("\n");
    io::stdout().flush().unwrap();

    self.simulation_time = start.elapsed();

    info!("Simulation completed in {:.2?} seconds.", self.simulation_time);

    if save {
      self.save_in_laamps_format("../output").unwrap();
    }
  }

  pub fn to_transfer_struct(&self) -> EngineDTO {
    EngineDTO {
      world: self.world.to_transfer_struct(),
      time_step: self.time_step,
      num_of_iterations: self.num_of_iterations,
    }
  }

  pub fn save_in_laamps_format(&self, path: &str) -> std::io::Result<()> {

    let start = Instant::now();
    let now: DateTime<Local> = Local::now();
    let time_string = now.format("%Y-%m-%d_%H-%M-%S").to_string();

    fs::create_dir(&format!("./{}/{}", path, time_string))?;

    let engine_dto = self.to_transfer_struct();

    for i in 0..engine_dto.num_of_iterations {
      if i % 100 == 0 {
        let mut result_string = String::new();
        result_string.push_str(&"ITEM: TIMESTEP\n".to_string());
        result_string.push_str(&format!("{}\n", i+1));

        result_string.push_str(&"ITEM: NUMBER OF ATOMS\n".to_string());
        let atom_container = engine_dto.world.atoms.get(i).unwrap();
        let num_of_atoms = atom_container.len();
        result_string.push_str(&format!("{}\n", num_of_atoms));

        result_string.push_str(&"ITEM: BOX BOUNDS pp pp pp\n".to_string());
        result_string.push_str(&format!("0.0 {}\n", engine_dto.world.box_x));
        result_string.push_str(&format!("0.0 {}\n", engine_dto.world.box_y));
        result_string.push_str(&format!("0.0 {}\n", engine_dto.world.box_z));

        result_string.push_str(&"ITEM: ATOMS id type x y z\n".to_string());

        for atom_dto in atom_container.iter() {
          result_string.push_str(&format!("{} {} {} {} {}\n", atom_dto.id, atom_dto.atom_type, atom_dto.x, atom_dto.y, atom_dto.z));
        }

        fs::write(&format!("./{}/{}/output_{}.dump", path, time_string, i+1), result_string)?;
      }
    }

    let mut kinetic_energy: Vec<f64> = Vec::with_capacity(engine_dto.num_of_iterations);
    // let mut potential_energy: Vec<f64> = Vec::with_capacity(engine_dto.num_of_iterations);
    let potential_energy = engine_dto.world.potential_energy;
    let mut total_energy: Vec<f64> = Vec::with_capacity(engine_dto.num_of_iterations);

    let mut forces: Vec<HashMap<usize, Vector3<f64>>> = Vec::new();

    for i in 0..engine_dto.num_of_iterations + 1 {
      let atom_container = engine_dto.world.atoms.get(i).unwrap();
      let mut current_forces: HashMap<usize, Vector3<f64>> = HashMap::new();
      let mut kinetic_energy_i = 0.;
      // let mut potential_energy_i = 0.;

      for atom_dto in atom_container.iter() {
        kinetic_energy_i += atom_dto.kinetic_energy;
        // potential_energy_i += atom_dto.potential_energy;

        current_forces.insert(atom_dto.id as usize, Vector3::new(atom_dto.force_x, atom_dto.force_y, atom_dto.force_z));
      }

      kinetic_energy.push(kinetic_energy_i);
      // potential_energy.push(potential_energy_i);
      let potential_energy_i = *potential_energy.get(i).unwrap();
      total_energy.push(kinetic_energy_i + potential_energy_i);
      forces.push(current_forces);
    }

    let mut wtr = Writer::from_path(&format!("./{}/{}/energy.csv", path, time_string))?;

    assert!(kinetic_energy.len() == potential_energy.len() && kinetic_energy.len() == total_energy.len());

    for i in 0..kinetic_energy.len() {
      wtr.write_record(&[
        format!("{}", i+1),
        format!("{}", kinetic_energy.get(i).unwrap()),
        format!("{}", potential_energy.get(i).unwrap()),
        format!("{}", total_energy.get(i).unwrap()),
      ])?;
    }

    wtr.flush()?;

    let mut wtr = Writer::from_path(&format!("./{}/{}/forces.csv", path, time_string))?;
    for i in 0..forces.len() {
      let force_map = forces.get(i).unwrap();

      for (id, force_vector) in force_map.iter() {
        wtr.write_record(&[
          format!("{}", i+1),
          format!("{}", id),
          format!("{}", force_vector.x),
          format!("{}", force_vector.y),
          format!("{}", force_vector.z),
        ])?;
      }
    }

    wtr.flush()?;

    let elapsed = start.elapsed();
    info!("Data saved in LAMMPS format in {:.2?} seconds.", elapsed);

    let mut wtr = Writer::from_path(&format!("./{}/{}/info.txt", path, time_string))?;
    wtr.write_record(&["Simulation Time", &format!("{:.2?} seconds", self.simulation_time)])?;
    wtr.write_record(&["Data Saving Time", &format!("{:.2?} seconds", elapsed)])?;
    wtr.flush()?;

    Ok(())
  }
}