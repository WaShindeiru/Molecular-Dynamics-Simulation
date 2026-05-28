use crate::data::types::AtomType;
use crate::particle::{Atom, CustomPathAtom, CustomVelocityAtom};
use crate::persistence::dto::atom::AtomDTO;
use crate::sim_core::world::computation::ForceComputationOperations;
use nalgebra::Vector3;
use std::any::Any;

#[derive(Debug, PartialEq, Eq, Clone, Copy)]
pub enum ParticleKind {
  Atom,
  CustomPathAtom,
  CustomVelocityAtom,
}

#[derive(Debug, PartialEq, Clone)]
pub enum Particle {
  Atom(Atom),
  CustomPathAtom(CustomPathAtom),
  CustomVelocityAtom(CustomVelocityAtom),
}

impl Particle {
  pub fn get_id(&self) -> usize {
    match self {
      Particle::Atom(atom) => atom.get_id(),
      Particle::CustomPathAtom(p) => p.get_id(),
      Particle::CustomVelocityAtom(p) => p.get_id(),
    }
  }

  pub fn get_type(&self) -> AtomType {
    match self {
      Particle::Atom(atom) => atom.get_type(),
      Particle::CustomPathAtom(p) => p.get_type(),
      Particle::CustomVelocityAtom(p) => p.get_type(),
    }
  }

  pub fn get_position(&self) -> &Vector3<f64> {
    match self {
      Particle::Atom(atom) => atom.get_position(),
      Particle::CustomPathAtom(p) => p.get_position(),
      Particle::CustomVelocityAtom(p) => p.get_position(),
    }
  }

  pub fn get_velocity(&self) -> &Vector3<f64> {
    match self {
      Particle::Atom(atom) => atom.get_velocity(),
      Particle::CustomPathAtom(p) => p.get_velocity(),
      Particle::CustomVelocityAtom(p) => p.get_velocity(),
    }
  }

  pub fn get_acceleration(&self) -> &Vector3<f64> {
    match self {
      Particle::Atom(atom) => atom.get_acceleration(),
      Particle::CustomPathAtom(p) => p.get_acceleration(),
      Particle::CustomVelocityAtom(p) => p.get_acceleration(),
    }
  }

  pub fn get_mass(&self) -> f64 {
    match self {
      Particle::Atom(atom) => atom.get_mass(),
      Particle::CustomPathAtom(p) => p.get_mass(),
      Particle::CustomVelocityAtom(p) => p.get_mass(),
    }
  }

  pub fn get_force(&self) -> &Vector3<f64> {
    match self {
      Particle::Atom(atom) => atom.get_force(),
      Particle::CustomPathAtom(p) => p.get_force(),
      Particle::CustomVelocityAtom(p) => p.get_force(),
    }
  }

  pub fn get_potential_energy(&self) -> f64 {
    match self {
      Particle::Atom(atom) => atom.get_potential_energy(),
      Particle::CustomPathAtom(p) => p.get_potential_energy(),
      Particle::CustomVelocityAtom(p) => p.get_potential_energy(),
    }
  }

  pub fn get_kinetic_energy(&self) -> f64 {
    match self {
      Particle::Atom(atom) => atom.get_kinetic_energy(),
      Particle::CustomPathAtom(p) => p.get_kinetic_energy(),
      Particle::CustomVelocityAtom(p) => p.get_kinetic_energy(),
    }
  }

  pub fn get_potential_gravity_energy(&self) -> f64 {
    match self {
      Particle::Atom(atom) => atom.get_potential_gravity_energy(),
      Particle::CustomPathAtom(p) => p.get_potential_gravity_energy(),
      Particle::CustomVelocityAtom(p) => p.get_potential_gravity_energy(),
    }
  }

  pub fn get_thermostat_work(&self) -> f64 {
    match self {
      Particle::Atom(atom) => atom.get_thermostat_work(),
      Particle::CustomPathAtom(p) => p.get_thermostat_work(),
      Particle::CustomVelocityAtom(p) => p.get_thermostat_work(),
    }
  }

  pub fn set_velocity(&mut self, velocity_: Vector3<f64>) {
    match self {
      Particle::Atom(atom) => atom.set_velocity(velocity_),
      Particle::CustomPathAtom(p) => p.set_velocity(velocity_),
      Particle::CustomVelocityAtom(p) => p.set_velocity(velocity_),
    }
  }

  pub fn set_potential_energy(&mut self, potential_energy_: f64) {
    match self {
      Particle::Atom(atom) => atom.set_potential_energy(potential_energy_),
      Particle::CustomPathAtom(p) => p.set_potential_energy(potential_energy_),
      Particle::CustomVelocityAtom(p) => p.set_potential_energy(potential_energy_),
    }
  }

  pub fn set_potential_gravity_energy(&mut self, potential_gravity_energy: f64) {
    match self {
      Particle::Atom(atom) => atom.set_potential_gravity_energy(potential_gravity_energy),
      Particle::CustomPathAtom(p) => p.set_potential_gravity_energy(potential_gravity_energy),
      Particle::CustomVelocityAtom(p) => p.set_potential_gravity_energy(potential_gravity_energy),
    }
  }

  pub fn set_force(&mut self, force_: Vector3<f64>) {
    match self {
      Particle::Atom(atom) => atom.set_force(force_),
      Particle::CustomPathAtom(p) => p.set_force(force_),
      Particle::CustomVelocityAtom(p) => p.set_force(force_),
    }
  }

  pub fn set_acceleration(&mut self, acceleration_: Vector3<f64>) {
    match self {
      Particle::Atom(atom) => atom.set_acceleration(acceleration_),
      Particle::CustomPathAtom(p) => p.set_acceleration(acceleration_),
      Particle::CustomVelocityAtom(p) => p.set_acceleration(acceleration_),
    }
  }

  pub fn set_thermostat_work(&mut self, thermostat_work: f64) {
    match self {
      Particle::Atom(atom) => atom.set_thermostat_work(thermostat_work),
      Particle::CustomPathAtom(p) => p.set_thermostat_work(thermostat_work),
      Particle::CustomVelocityAtom(p) => p.set_thermostat_work(thermostat_work),
    }
  }

  pub fn get_iteration(&self) -> usize {
    match self {
      Particle::Atom(atom) => atom.get_iteration(),
      Particle::CustomPathAtom(p) => p.get_iteration(),
      Particle::CustomVelocityAtom(p) => p.get_iteration(),
    }
  }

  pub fn set_iteration(&mut self, iteration_: usize) {
    match self {
      Particle::Atom(atom) => atom.set_iteration(iteration_),
      Particle::CustomPathAtom(p) => p.set_iteration(iteration_),
      Particle::CustomVelocityAtom(p) => p.set_iteration(iteration_),
    }
  }

  pub fn update_position(&mut self, position_: Vector3<f64>) {
    match self {
      Particle::Atom(atom) => atom.update_position(position_),
      Particle::CustomPathAtom(p) => p.update_position(position_),
      Particle::CustomVelocityAtom(p) => p.update_position(position_),
    }
  }

  pub fn to_transfer_struct(&self) -> AtomDTO {
    match self {
      Particle::Atom(atom) => atom.to_transfer_struct(),
      Particle::CustomPathAtom(p) => p.to_transfer_struct(),
      Particle::CustomVelocityAtom(p) => p.to_transfer_struct(),
    }
  }

  pub fn custom_clone(&self) -> Particle {
    match self {
      Particle::Atom(atom) => Particle::Atom(atom.clone()),
      Particle::CustomPathAtom(p) => Particle::CustomPathAtom(p.clone()),
      Particle::CustomVelocityAtom(p) => Particle::CustomVelocityAtom(p.clone()),
    }
  }

  pub fn reset_clone(&self) -> Particle {
    match self {
      Particle::Atom(atom) => Particle::Atom(atom.reset_clone()),
      Particle::CustomPathAtom(p) => Particle::CustomPathAtom(p.reset_clone()),
      Particle::CustomVelocityAtom(p) => Particle::CustomVelocityAtom(p.reset_clone()),
    }
  }
}

impl AsRef<Particle> for Particle {
  fn as_ref(&self) -> &Particle {
    self
  }
}

impl ForceComputationOperations for Particle {
  fn as_any(&self) -> &dyn Any {
    self
  }

  fn get_id(&self) -> usize {
    Particle::get_id(self)
  }

  fn get_position(&self) -> Vector3<f64> {
    *Particle::get_position(self)
  }

  fn get_type(&self) -> AtomType {
    Particle::get_type(self)
  }

  fn get_mass(&self) -> f64 {
    Particle::get_mass(self)
  }

  fn prototype_clone(&self) -> Box<dyn ForceComputationOperations> {
    Box::new(self.clone())
  }
}

impl Particle {
  pub fn is_custom_velocity_atom(&self) -> bool {
    matches!(self, Particle::CustomVelocityAtom(_))
  }
}

pub fn compute_kinetic_energy(particles: &Vec<Particle>) -> f64 {
  let mut kinetic_energy = 0.;

  for particle in particles {
    kinetic_energy += particle.get_mass() * particle.get_velocity().magnitude().powi(2) / 2.0;
  }

  kinetic_energy
}
