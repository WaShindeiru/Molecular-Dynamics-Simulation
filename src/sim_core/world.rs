pub mod verlet;
pub mod euler;
pub mod verlet_nose_hoover;
pub mod persist;
pub mod integration;

use log::info;
use nalgebra::Vector3;
use crate::particle::{Particle, SimpleAtomContainer};

use integration::IntegrationAlgorithm;

pub struct World {
  atoms: Vec<SimpleAtomContainer>,
  atom_count: usize,
  size: Vector3<f64>,
  current_iteration: usize,
  current_index: usize,
  thermostat_epsilon: Vec<f64>,
  thermostat_work_total: f64,

  max_iteration_till_reset: usize,
  reset_counter: usize,
  number_of_resets: usize,

  frame_iteration_count: usize,
  frame_iteration_count_current_iteration: usize,

  save: bool,
  save_path: String,
  save_laamps: bool,
  save_verbose: bool,
}

fn get_index_for_iteration(current_iteration: usize, max_iteration_till_reset: usize, reset_count: usize) -> usize {
  let mut index = current_iteration;
  if reset_count > 0 {
    index = index - (max_iteration_till_reset - 1);
  }
  index = index - (max_iteration_till_reset) * (reset_count - 1);
  index
}

impl World {
  // pub fn new(size: Vector3<f64>) -> Self {
  //   World {
  //     atoms: Vec::new(),
  //     atom_count: 0,
  //     size,
  //     current_iteration: 0,
  //     thermostat_epsilon: vec![0.],
  //   }
  // }

  pub fn new_from_atoms(atoms: Vec<Particle>, size: Vector3<f64>, max_iteration_till_reset: usize,
                        frame_iteration_count: usize, save: bool, save_path: String,
                        save_laamps: bool, save_verbose: bool) -> Self {
    let atom_count = atoms.len();

    let atom_container = SimpleAtomContainer::new_from_atoms(atoms);
    let mut atoms: Vec<SimpleAtomContainer> = Vec::with_capacity(max_iteration_till_reset);
    atoms.push(atom_container);

    let mut thermostat_epsilon: Vec<f64> = Vec::with_capacity(max_iteration_till_reset);
    thermostat_epsilon.push(0.);

    World {
      atoms,
      atom_count,
      size,
      current_iteration: 0,
      current_index: 0,
      thermostat_epsilon,
      thermostat_work_total: 0.,

      max_iteration_till_reset,
      reset_counter: 1,
      number_of_resets: 0,

      frame_iteration_count,
      frame_iteration_count_current_iteration: 0,

      save,
      save_path,
      save_laamps,
      save_verbose,
    }
  }

  pub fn update(&mut self, params: &IntegrationAlgorithm, time_step: f64, next_iteration: usize) {
    assert!(self.reset_counter <= self.max_iteration_till_reset);
    assert_eq!(self.atoms.len() - 1, self.current_index);

    if self.reset_counter == self.max_iteration_till_reset {
      let use_thermostat: bool;
      match params {
        IntegrationAlgorithm::NoseHooverVerlet {..} => {
          use_thermostat = true;
        }
        _ => {
          use_thermostat = false;
        }
      }

      if self.save {
        self.save(use_thermostat).unwrap();
      }
      self.reset_world();
    }

    match params {
      IntegrationAlgorithm::SemiImplicitEuler => {
        self.update_semi_implicit_euler(time_step, next_iteration);
      }
      IntegrationAlgorithm::VelocityVerlet => {
        self.update_verlet(time_step, next_iteration);
      }
      IntegrationAlgorithm::NoseHooverVerlet {
        desired_temperature,
        q_effective_mass } => {
        self.update_verlet_nose_hoover(time_step, next_iteration, *desired_temperature, *q_effective_mass);
      }
    }

    self.current_index += 1;
    self.reset_counter += 1;
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
}
