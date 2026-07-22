use nalgebra::base::Vector3;
use rand::distr::Uniform;
use rand_distr::Distribution;
use std::sync::{Mutex, OnceLock};

use crate::data::constants::{ATOMIC_MASS_C, ATOMIC_MASS_FE};
use crate::data::types::AtomType;
use crate::particle::Particle;
use crate::particle::custom_path_atom::CustomPathAtom;
use crate::particle::custom_velocity_atom::CustomVelocityAtom;
use crate::particle::particle::ParticleKind;
use crate::persistence::dto::atom::AtomDTO;

#[derive(Debug, PartialEq, Clone)]
pub struct Atom {
  id: usize,
  iteration: usize,
  type_: AtomType,
  mass: f64,

  position: Vector3<f64>,
  velocity: Vector3<f64>,

  acceleration: Vector3<f64>,
  force: Vector3<f64>,

  kinetic_energy: f64,
  potential_energy: f64,
  thermostat_work: f64,

  potential_gravity_energy: f64,
}

impl Atom {
  pub fn get_id(&self) -> usize {
    self.id
  }

  pub fn get_iteration(&self) -> usize {
    self.iteration
  }

  pub fn get_type(&self) -> AtomType {
    self.type_
  }

  pub fn get_position(&self) -> &Vector3<f64> {
    &self.position
  }

  pub fn get_velocity(&self) -> &Vector3<f64> {
    &self.velocity
  }

  pub fn get_acceleration(&self) -> &Vector3<f64> {
    &self.acceleration
  }

  pub fn get_mass(&self) -> f64 {
    self.mass
  }

  pub fn get_force(&self) -> &Vector3<f64> {
    &self.force
  }

  pub fn get_potential_energy(&self) -> f64 {
    self.potential_energy
  }

  pub fn get_kinetic_energy(&self) -> f64 {
    self.kinetic_energy
  }

  pub fn get_potential_gravity_energy(&self) -> f64 {
    self.potential_gravity_energy
  }

  pub fn get_thermostat_work(&self) -> f64 {
    self.thermostat_work
  }

  pub fn set_velocity(&mut self, velocity_: Vector3<f64>) {
    self.velocity = velocity_;
    self.kinetic_energy = self.mass * velocity_.magnitude_squared() / 2.0;
  }

  pub fn set_potential_energy(&mut self, potential_energy_: f64) {
    self.potential_energy = potential_energy_;
  }

  pub fn set_potential_gravity_energy(&mut self, potential_gravity_energy: f64) {
    self.potential_gravity_energy = potential_gravity_energy;
  }

  pub fn set_force(&mut self, force_: Vector3<f64>) {
    self.force = force_;
  }

  pub fn set_acceleration(&mut self, acceleration_: Vector3<f64>) {
    self.acceleration = acceleration_;
  }

  pub fn set_thermostat_work(&mut self, thermostat_work: f64) {
    self.thermostat_work = thermostat_work
  }

  pub fn set_iteration(&mut self, iteration_: usize) {
    self.iteration = iteration_;
  }

  pub fn update_position(&mut self, position_: Vector3<f64>) {
    self.position = position_;
  }

  pub fn set_position(&mut self, position_: Vector3<f64>) {
    self.position = position_;
  }

  pub fn reset_clone(&self) -> Atom {
    Atom {
      id: self.id,
      iteration: self.iteration,
      type_: self.type_,
      mass: self.mass,

      position: Vector3::new(0., 0., 0.),
      velocity: self.velocity, // TODO: reconsider this

      acceleration: Vector3::new(0., 0., 0.),
      force: Vector3::new(0., 0., 0.),

      kinetic_energy: self.kinetic_energy,
      potential_energy: 0.,
      thermostat_work: 0.,
      potential_gravity_energy: 0.,
    }
  }

  pub fn to_transfer_struct(&self) -> AtomDTO {
    AtomDTO {
      id: self.id,
      iteration: self.iteration,
      kind: ParticleKind::Atom,
      atom_type: self.type_,
      position: self.position,
      velocity: self.velocity,
      velocity_manager_id: None,
      control_velocity_manager_id: None,
      kinetic_energy: self.kinetic_energy,
      potential_energy: self.potential_energy,
      thermostat_work: self.thermostat_work,
      potential_gravity_energy: self.potential_gravity_energy,
      p_control_energy: 0.0,
      force: self.force,
    }
  }
}

impl Atom {
  pub fn new(
    id: usize,
    type_: AtomType,
    mass: f64,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
    force: Vector3<f64>,
    acceleration: Vector3<f64>,
    potential_energy: f64,
  ) -> Self {
    Atom {
      id,
      iteration: 0,
      type_,
      mass,
      position,
      velocity,
      force,
      acceleration,
      kinetic_energy: mass * velocity.magnitude_squared() / 2.0,
      potential_energy,
      potential_gravity_energy: 0.,
      thermostat_work: 0.0,
    }
  }

  pub fn new_custom_iteration(
    id: usize,
    iteration: usize,
    type_: AtomType,
    mass: f64,
    position: Vector3<f64>,
    velocity: Vector3<f64>,
    force: Vector3<f64>,
    acceleration: Vector3<f64>,
    potential_energy: f64,
  ) -> Self {
    Atom {
      id,
      iteration,
      type_,
      mass,
      position,
      velocity,
      force,
      acceleration,
      kinetic_energy: mass * velocity.magnitude_squared() / 2.0,
      potential_energy,
      potential_gravity_energy: 0.,
      thermostat_work: 0.0,
    }
  }
}

struct AtomFactory {
  counter: usize,
  potential_gravity_max: f64,
  z_max: f64,
}

impl AtomFactory {
  fn new(potential_gravity_max: f64, z_max: f64) -> Self {
    AtomFactory {
      counter: 0usize,
      potential_gravity_max,
      z_max,
    }
  }

  fn get_atom(
    &mut self,
    atom: AtomType,
    position_: Vector3<f64>,
    velocity_: Vector3<f64>,
  ) -> Particle {
    let result = match atom {
      AtomType::C | AtomType::C_nanotube | AtomType::C_nanotube_static => Atom {
        id: self.counter,
        iteration: 0,
        type_: atom,
        mass: ATOMIC_MASS_C,
        position: position_,
        velocity: velocity_,
        force: Vector3::new(0.0, 0.0, 0.0),
        acceleration: Vector3::new(0.0, 0.0, 0.0),
        kinetic_energy: ATOMIC_MASS_C * velocity_.magnitude_squared() / 2.0,
        potential_energy: 0.0,
        thermostat_work: 0.0,
        potential_gravity_energy: self.potential_gravity_max * position_.z * ATOMIC_MASS_C
          / self.z_max,
      },
      AtomType::Fe => Atom {
        type_: AtomType::Fe,
        iteration: 0,
        id: self.counter,
        mass: ATOMIC_MASS_FE,
        position: position_,
        velocity: velocity_,
        force: Vector3::new(0.0, 0.0, 0.0),
        acceleration: Vector3::new(0.0, 0.0, 0.0),
        kinetic_energy: ATOMIC_MASS_FE * velocity_.magnitude_squared() / 2.0,
        potential_energy: 0.0,
        thermostat_work: 0.0,
        potential_gravity_energy: self.potential_gravity_max * position_.z * ATOMIC_MASS_FE
          / self.z_max,
      },
    };
    self.counter = self.counter + 1;

    Particle::Atom(result)
  }

  fn get_atom_custom_path(&mut self, atom: AtomType, path: Vec<Vector3<f64>>) -> Particle {
    let result = match atom {
      AtomType::C | AtomType::C_nanotube | AtomType::C_nanotube_static => CustomPathAtom::new(
        self.counter,
        atom,
        ATOMIC_MASS_C,
        path,
        self.potential_gravity_max,
        self.z_max,
      ),
      AtomType::Fe => CustomPathAtom::new(
        self.counter,
        AtomType::Fe,
        ATOMIC_MASS_FE,
        path,
        self.potential_gravity_max,
        self.z_max,
      ),
    };
    self.counter = self.counter + 1;

    Particle::CustomPathAtom(result)
  }

  fn get_atom_custom_velocity(
    &mut self,
    atom: AtomType,
    position: Vector3<f64>,
    particle_velocity_manager_id: usize,
  ) -> Particle {
    let result = match atom {
      AtomType::C | AtomType::C_nanotube | AtomType::C_nanotube_static => {
        CustomVelocityAtom::new(self.counter, atom, ATOMIC_MASS_C, position, particle_velocity_manager_id)
      }
      AtomType::Fe => CustomVelocityAtom::new(
        self.counter,
        AtomType::Fe,
        ATOMIC_MASS_FE,
        position,
        particle_velocity_manager_id,
      ),
    };
    self.counter = self.counter + 1;

    Particle::CustomVelocityAtom(result)
  }
}

pub struct SafeAtomFactory {
  potential_gravity_max: f64,
  z_max: f64,
  inner: Mutex<AtomFactory>,
}

static SAFE_ATOM_FACTORY_INSTANCE: OnceLock<SafeAtomFactory> = OnceLock::new();

impl SafeAtomFactory {
  pub fn new(potential_gravity_max: f64, z_max: f64) -> &'static Self {
    let instance = SAFE_ATOM_FACTORY_INSTANCE.get_or_init(|| Self {
      potential_gravity_max,
      z_max,
      inner: Mutex::new(AtomFactory::new(potential_gravity_max, z_max)),
    });

    assert!(
      (instance.potential_gravity_max - potential_gravity_max).abs() <= f64::EPSILON
        && (instance.z_max - z_max).abs() <= f64::EPSILON,
      "SafeAtomFactory singleton already initialized with potential_gravity_max={} and z_max={}, but requested potential_gravity_max={} and z_max={}",
      instance.potential_gravity_max,
      instance.z_max,
      potential_gravity_max,
      z_max
    );

    instance
  }

  pub fn reset_for_testing() {
    // SAFETY: Only called from tests that hold the FACTORY_TEST_LOCK, which serialises all
    // callers. No thread can be concurrently reading SAFE_ATOM_FACTORY_INSTANCE while this
    // runs because they are all blocked on the same lock.
    unsafe {
      let ptr = &raw const SAFE_ATOM_FACTORY_INSTANCE as *mut OnceLock<SafeAtomFactory>;
      std::mem::replace(&mut *ptr, OnceLock::new());
    }
  }

  pub fn get_atom(
    &self,
    atom: AtomType,
    position_: Vector3<f64>,
    velocity_: Vector3<f64>,
  ) -> Particle {
    let mut factory = self.inner.lock().unwrap();
    factory.get_atom(atom, position_, velocity_)
  }

  pub fn get_atom_custom_path(&self, atom: AtomType, path: Vec<Vector3<f64>>) -> Particle {
    let mut factory = self.inner.lock().unwrap();
    factory.get_atom_custom_path(atom, path)
  }

  pub fn get_atom_custom_velocity(
    &self,
    atom: AtomType,
    position: Vector3<f64>,
    particle_velocity_manager_id: usize,
  ) -> Particle {
    let mut factory = self.inner.lock().unwrap();
    factory.get_atom_custom_velocity(atom, position, particle_velocity_manager_id)
  }

  pub fn get_atom_random(&self, atom: AtomType, lower_bound: f64, upper_bound: f64) -> Particle {
    let mut factory = self.inner.lock().unwrap();

    let range = Uniform::new(lower_bound, upper_bound).unwrap();
    let range_vel = Uniform::new(-1., 1.).unwrap();
    let mut rng = rand::rng();

    let position_ = Vector3::new(
      range.sample(&mut rng),
      range.sample(&mut rng),
      range.sample(&mut rng),
    );
    let velocity_ = Vector3::new(
      range_vel.sample(&mut rng),
      range_vel.sample(&mut rng),
      range_vel.sample(&mut rng),
    );
    factory.get_atom(atom, position_, velocity_)
  }
}
