use std::cmp::Reverse;
use std::collections::BinaryHeap;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::data::constants::{ATOMIC_MASS_C, ATOMIC_MASS_FE};
use crate::data::types::AtomType;
use ParticleKind::*;

use crate::particle::particle::{Particle, ParticleKind};
use crate::persistence::dto::atom::AtomDTO;

#[derive(PartialEq)]
pub struct VelocityParticle {
  pub id: usize,
  pub iteration: usize,
  pub kind: ParticleKind,
  pub type_: AtomType,
  pub mass: f64,
  pub position: Vector3<f64>,
  pub velocity: Vector3<f64>,
  pub kinetic_energy: f64,
}

impl Eq for VelocityParticle {}

impl PartialOrd for VelocityParticle {
  fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
    Some(self.cmp(other))
  }
}

impl Ord for VelocityParticle {
  fn cmp(&self, other: &Self) -> std::cmp::Ordering {
    match (&self.kind, &other.kind) {
      (Atom, Atom) => self.velocity.magnitude().total_cmp(&other.velocity.magnitude()),
    
      (Atom, CustomPathAtom) | (Atom, CustomVelocityAtom) => std::cmp::Ordering::Greater,
      (CustomPathAtom, Atom) | (CustomVelocityAtom, Atom) => std::cmp::Ordering::Less,

      (CustomPathAtom, CustomPathAtom)
      | (CustomVelocityAtom, CustomVelocityAtom)
      | (CustomPathAtom, CustomVelocityAtom)
      | (CustomVelocityAtom, CustomPathAtom) => std::cmp::Ordering::Equal,
    }
  }
}

impl From<&Particle> for VelocityParticle {
  fn from(p: &Particle) -> Self {
    let kind = match p {
      Particle::Atom(_) => ParticleKind::Atom,
      Particle::CustomPathAtom(_) => ParticleKind::CustomPathAtom,
      Particle::CustomVelocityAtom(_) => ParticleKind::CustomVelocityAtom,
    };
    
    VelocityParticle {
      id: p.get_id(),
      iteration: p.get_iteration(),
      kind,
      type_: p.get_type(),
      mass: p.get_mass(),
      position: *p.get_position(),
      velocity: *p.get_velocity(),
      kinetic_energy: p.get_kinetic_energy(),
    }
  }
}

impl From<Particle> for VelocityParticle {
  fn from(p: Particle) -> Self {
    VelocityParticle::from(&p)
  }
}

impl From<AtomDTO> for VelocityParticle {
  fn from(dto: AtomDTO) -> Self {
    let mass = match dto.atom_type {
      AtomType::C | AtomType::C_nanotube => ATOMIC_MASS_C,
      AtomType::Fe => ATOMIC_MASS_FE,
    };
    let kinetic_energy = 0.5 * mass * dto.velocity.magnitude_squared();
    VelocityParticle {
      id: dto.id,
      iteration: dto.iteration,
      kind: dto.kind,
      type_: dto.atom_type,
      mass,
      position: dto.position,
      velocity: dto.velocity,
      kinetic_energy,
    }
  }
}

impl From<Arc<Particle>> for VelocityParticle {
  fn from(p: Arc<Particle>) -> Self {
    VelocityParticle::from(p.as_ref())
  }
}

pub struct VelocityHeap {
  heap: BinaryHeap<Reverse<VelocityParticle>>,
  capacity: usize,
}

impl VelocityHeap {
  pub fn new(capacity: usize) -> VelocityHeap {
    VelocityHeap {
      heap: BinaryHeap::new(),
      capacity,
    }
  }

  pub fn try_insert<T: Into<VelocityParticle>>(&mut self, value: T) {
    let vp = value.into();
    if self.heap.len() < self.capacity {
      self.heap.push(Reverse(vp));

    } else if let Some(Reverse(min)) = self.heap.peek() {
      if vp > *min {
        self.heap.pop();
        self.heap.push(Reverse(vp));
      }
    }
  }

  pub fn drain_descending(&mut self) -> impl Iterator<Item = VelocityParticle> {
    let mut items: Vec<_> = self.heap.drain().map(|Reverse(vp)| vp).collect();
    items.sort_unstable_by(|a, b| b.cmp(a));
    items.into_iter()
  }
}