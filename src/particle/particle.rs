use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::output::AtomDTO;
use crate::particle::{Atom, CustomPathAtom};

pub enum Particle {
  Atom(Atom),
  CustomPathAtom(CustomPathAtom)
}

impl Particle {
  pub fn get_id(&self) -> u64 {
    match self {
      Particle::Atom(atom) => atom.get_id(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_id(),
    }
  }

  pub fn get_type(&self) -> &AtomType {
    match self {
      Particle::Atom(atom) => atom.get_type(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_type(),
    }
  }

  pub fn get_position(&self) -> &Vector3<f64> {
    match self {
      Particle::Atom(atom) => atom.get_position(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_position(),
    }
  }

  pub fn get_velocity(&self) -> &Vector3<f64> {
    match self {
      Particle::Atom(atom) => atom.get_velocity(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_velocity(),
    }
  }

  pub fn get_acceleration(&self) -> &Vector3<f64> {
    match self {
      Particle::Atom(atom) => atom.get_acceleration(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_acceleration(),
    }
  }
  
  pub fn get_mass(&self) -> f64 {
    match self {
      Particle::Atom(atom) => atom.get_mass(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_mass(),
    }
  }
  
  pub fn get_force(&self) -> &Vector3<f64> {
    match self {
      Particle::Atom(atom) => atom.get_force(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_force(),
    }
  }
  
  pub fn get_potential_energy(&self) -> f64 {
    match self {
      Particle::Atom(atom) => atom.get_potential_energy(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_potential_energy(),
    }
  }

  pub fn get_thermostat_work(&self) -> f64 {
    match self {
      Particle::Atom(atom) => atom.get_thermostat_work(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.get_thermostat_work(),
    }
  }
  
  pub fn set_velocity(&mut self, velocity_: Vector3<f64>) {
    match self {
      Particle::Atom(atom) => atom.set_velocity(velocity_),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.set_velocity(velocity_),
    }
  }
  
  pub fn set_potential_energy(&mut self, potential_energy_: f64) {
    match self {
      Particle::Atom(atom) => atom.set_potential_energy(potential_energy_),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.set_potential_energy(potential_energy_),
    }
  }
  
  pub fn set_force(&mut self, force_: Vector3<f64>) {
    match self {
      Particle::Atom(atom) => atom.set_force(force_),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.set_force(force_),
    }
  }
  
  pub fn set_acceleration(&mut self, acceleration_: Vector3<f64>) {
    match self {
      Particle::Atom(atom) => atom.set_acceleration(acceleration_),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.set_acceleration(acceleration_),
    }
  }

  pub fn set_thermostat_work(&mut self, thermostat_work: f64) {
    match self {
      Particle::Atom(atom) => atom.set_thermostat_work(thermostat_work),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.set_thermostat_work(thermostat_work),
    }
  }

  pub fn set_iteration(&mut self, iteration_: usize) {
    match self {
      Particle::Atom(atom) => atom.set_iteration(iteration_),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.set_iteration(iteration_),
    }
  }
  
  pub fn update_position(&mut self, position_: Vector3<f64>) {
    match self {
      Particle::Atom(atom) => atom.update_position(position_),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.update_position(position_),
    }
  }
  
  pub fn to_transfer_struct(&self) -> AtomDTO {
    match self {
      Particle::Atom(atom) => atom.to_transfer_struct(),
      Particle::CustomPathAtom(custom_path_atom) => custom_path_atom.to_transfer_struct(),
    }
  }
  
  pub fn custom_clone(&self) -> Particle {
    match self {
      Particle::Atom(atom) => Particle::Atom(atom.clone()),
      Particle::CustomPathAtom(custom_path_atom) => Particle::CustomPathAtom(custom_path_atom.clone()),
    }
  }
}

pub fn compute_kinetic_energy(particles: &Vec<Particle>) -> f64 {
  let mut kinetic_energy = 0.;

  for particle in particles {
    kinetic_energy += particle.get_mass() * particle.get_velocity().magnitude().powi(2) / 2.0;
  }

  kinetic_energy
}