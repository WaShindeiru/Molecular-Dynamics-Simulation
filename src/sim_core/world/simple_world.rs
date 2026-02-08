use std::io;
use log::info;
use nalgebra::Vector3;
use crate::output::{AtomDTO, WorldDTO};
use crate::sim_core::world::saver::PartialWorldSaver;
use crate::particle::{Particle, SimpleAtomContainer};
use crate::sim_core::world::get_index_for_iteration;
use crate::sim_core::world::integration::{IntegrationAlgorithm, IntegrationAlgorithmParams};
use crate::sim_core::world::saver::SaveOptions;

pub mod integration;

pub struct SimpleWorld {
  atoms: Vec<SimpleAtomContainer>,
  atom_count: usize,
  size: Vector3<f64>,
  current_iteration: usize,
  current_index: usize,
  integration_algorithm: IntegrationAlgorithm,

  thermostat_epsilon: Vec<f64>,
  thermostat_work_total: f64,

  max_iteration_till_reset: usize,
  reset_counter: usize,
  number_of_resets: usize,

  frame_iteration_count: usize,

  save_options: SaveOptions,
  world_saver: PartialWorldSaver,
}

impl SimpleWorld {
  pub fn new_from_atoms(atoms: Vec<Particle>, size: Vector3<f64>, max_iteration_till_reset: usize,
                        frame_iteration_count: usize, integration_algorithm: IntegrationAlgorithm, save_options: SaveOptions) -> Self {
    let atom_count = atoms.len();

    let atom_container = SimpleAtomContainer::new_from_atoms(atoms);
    let mut atoms: Vec<SimpleAtomContainer> = Vec::with_capacity(max_iteration_till_reset);
    atoms.push(atom_container);

    let mut thermostat_epsilon: Vec<f64> = Vec::with_capacity(max_iteration_till_reset);
    thermostat_epsilon.push(0.);

    SimpleWorld {
      atoms,
      atom_count,
      size,
      current_iteration: 0,
      current_index: 0,
      integration_algorithm,

      thermostat_epsilon,
      thermostat_work_total: 0.,

      max_iteration_till_reset,
      reset_counter: 1,
      number_of_resets: 0,

      frame_iteration_count,

      save_options: save_options.clone(),
      world_saver: PartialWorldSaver::new(save_options),
    }
  }

  pub fn save(&mut self) -> io::Result<()> {
    if self.save_options.save {
      let world = self.to_transfer_struct();
      self.world_saver.persist(&world)?;
    }
    Ok(())
  }

  pub fn update(&mut self, params: &IntegrationAlgorithmParams, time_step: f64, next_iteration: usize) {
    assert!(self.reset_counter <= self.max_iteration_till_reset);
    assert_eq!(self.atoms.len() - 1, self.current_index);

    if self.reset_counter == self.max_iteration_till_reset {
      self.save().unwrap();
      self.reset_world();
    }

    match params {
      IntegrationAlgorithmParams::SemiImplicitEuler => {
        self.update_semi_implicit_euler(time_step, next_iteration);
      }
      IntegrationAlgorithmParams::VelocityVerlet => {
        self.update_verlet(time_step, next_iteration);
      }
      IntegrationAlgorithmParams::NoseHooverVerlet {
        desired_temperature,
        q_effective_mass } => {
        self.update_verlet_nose_hoover(time_step, next_iteration, *desired_temperature, *q_effective_mass);
      }
    }

    self.current_index += 1;
    self.reset_counter += 1;
  }

  pub fn get_size(&self) -> &Vector3<f64> {
    &self.size
  }

  pub fn reset_world(&mut self) {
    info!("Resetting world: {}", self.number_of_resets);

    let new_number_of_resets = self.number_of_resets + 1;
    let new_index = get_index_for_iteration(self.current_iteration, self.max_iteration_till_reset, new_number_of_resets);

    let mut new_atoms: Vec<SimpleAtomContainer> = Vec::with_capacity(self.max_iteration_till_reset + 1);
    if let Some(last_atom_container) = self.atoms.pop() {
      new_atoms.push(last_atom_container);
    } else {
      panic!("world is empty!");
    }

    let mut new_thermostat_epsilon: Vec<f64> = Vec::with_capacity(self.max_iteration_till_reset + 1);
    if let Some(last_epsilon) = self.thermostat_epsilon.pop() {
      new_thermostat_epsilon.push(last_epsilon);
    } else {
      panic!("world is empty!");
    }

    self.atoms = new_atoms;
    self.current_index = new_index;
    self.thermostat_epsilon = new_thermostat_epsilon;
    self.number_of_resets = new_number_of_resets;
    self.reset_counter = 0;
  }

  pub fn apply_boundary_constraint(&self, mut atom: Particle) -> Particle {
    atom = match atom.get_position().x {
      x if x < 0.0 => {
        atom.set_velocity(Vector3::new(-atom.get_velocity().x, atom.get_velocity().y, atom.get_velocity().z));
        atom.update_position(Vector3::new(-atom.get_position().x, atom.get_position().y, atom.get_position().z));
        atom
      }
      x if x > self.size.x => {
        atom.set_velocity(Vector3::new(-atom.get_velocity().x, atom.get_velocity().y, atom.get_velocity().z));
        atom.update_position(Vector3::new( 2. * self.size.x - atom.get_position().x, atom.get_position().y, atom.get_position().z));
        atom
      },
      _ => atom
    };

    atom = match atom.get_position().y {
      y if y < 0.0 => {
        atom.set_velocity(Vector3::new(atom.get_velocity().x, -atom.get_velocity().y, atom.get_velocity().z));
        atom.update_position(Vector3::new(atom.get_position().x, -atom.get_position().y, atom.get_position().z));
        atom
      }
      y if y > self.size.y => {
        atom.set_velocity(Vector3::new(atom.get_velocity().x, -atom.get_velocity().y, atom.get_velocity().z));
        atom.update_position(Vector3::new(atom.get_position().x, 2. * self.size.y - atom.get_position().y, atom.get_position().z));
        atom
      },
      _ => atom
    };

    atom = match atom.get_position().z {
      z if z < 0.0 => {
        atom.set_velocity(Vector3::new(atom.get_velocity().x, atom.get_velocity().y, -atom.get_velocity().z));
        atom.update_position(Vector3::new(atom.get_position().x, atom.get_position().y, -atom.get_position().z));
        atom
      }
      z if z > self.size.z => {
        atom.set_velocity(Vector3::new(atom.get_velocity().x, atom.get_velocity().y, -atom.get_velocity().z));
        atom.update_position(Vector3::new(atom.get_position().x, atom.get_position().y, 2. * self.size.z - atom.get_position().z));
        atom
      },
      _ => atom
    };

    atom
  }

  pub fn to_transfer_struct(&self) -> WorldDTO {
    let mut all_atoms_dto: Vec<Vec<AtomDTO>> = Vec::with_capacity(self.atoms.len());
    let mut potential_energies: Vec<f64> = Vec::with_capacity(self.atoms.len());

    let lower_index: usize;

    if self.number_of_resets > 0 {
      lower_index = 1;
    } else {
      lower_index = 0;
    }

    for iteration in lower_index..self.atoms.len() {
      let atom_container = self.atoms.get(iteration).unwrap();
      let mut atoms_dto: Vec<AtomDTO> = Vec::with_capacity(atom_container.len());

      for atom in atom_container.get_atoms().iter() {
        atoms_dto.push(atom.to_transfer_struct());
      }

      all_atoms_dto.push(atoms_dto);
      potential_energies.push(atom_container.get_potential_energy());
    }

    WorldDTO {
      num_of_atoms: self.atom_count,
      atoms: all_atoms_dto,
      potential_energy: potential_energies,
      thermostat_epsilon: self.thermostat_epsilon.clone(),
      box_x: self.size.x,
      box_y: self.size.y,
      box_z: self.size.z,
      integration_algorithm: self.integration_algorithm.clone(),

      num_of_world_iterations: self.reset_counter,
      number_of_resets: self.number_of_resets,
      max_iteration_till_reset: self.max_iteration_till_reset,

      frame_iteration_count: self.frame_iteration_count,
    }
  }
}