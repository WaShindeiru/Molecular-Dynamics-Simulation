pub mod verlet;
pub mod euler;
pub mod verlet_nose_hoover;

use nalgebra::Vector3;
use crate::output::{AtomDTO, WorldDTO};
use crate::particle::{Particle, SimpleAtomContainer};

pub struct World {
  atoms: Vec<SimpleAtomContainer>,
  atom_count: usize,
  size: Vector3<f64>,
  current_iteration: usize,
  thermostat_epsilon: Vec<f64>,
}

impl World {
  pub fn new(size: Vector3<f64>) -> Self {
    World {
      atoms: Vec::new(),
      atom_count: 0,
      size,
      current_iteration: 0,
      thermostat_epsilon: vec![0.],
    }
  }

  pub fn new_from_atoms(atoms: Vec<Particle>, size: Vector3<f64>) -> Self {
    let atom_count = atoms.len();
    let atom_container = SimpleAtomContainer::new_from_atoms(atoms);
    World {
      atoms: vec![atom_container],
      atom_count,
      size,
      current_iteration: 0,
      thermostat_epsilon: vec![0.],
    }
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

    for iteration in 0..self.atoms.len() {
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
    }
  }
}
