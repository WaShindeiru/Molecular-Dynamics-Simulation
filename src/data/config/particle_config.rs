use crate::data::types::AtomType;
use crate::particle::Particle;
use nalgebra::Vector3;

#[derive(Clone)]
pub struct VelocityScheduleConfig {
  pub particle_id: usize,
  /// Iteration-velocity pairs; must include an entry for iteration 0.
  pub velocities: Vec<(usize, Vector3<f64>)>,
}

#[derive(Clone)]
pub struct ParticleConfig {
  pub atoms: Vec<Particle>,
  pub num_of_atoms: usize,
  pub num_of_carbon_atoms: usize,
  pub num_of_iron_atoms: usize,
  pub velocity_schedules: Vec<VelocityScheduleConfig>,
}

impl ParticleConfig {
  pub fn new(atoms: Vec<Particle>) -> Self {
    Self::new_with_schedules(atoms, vec![])
  }

  pub fn new_with_schedules(
    atoms: Vec<Particle>,
    velocity_schedules: Vec<VelocityScheduleConfig>,
  ) -> Self {
    let (num_of_carbon_atoms, num_of_iron_atoms) =
      atoms
        .iter()
        .fold((0usize, 0usize), |acc, atom| match atom.get_type() {
          AtomType::C | AtomType::C_nanotube => (acc.0 + 1, acc.1),
          AtomType::Fe => (acc.0, acc.1 + 1),
        });

    ParticleConfig {
      num_of_atoms: atoms.len(),
      num_of_carbon_atoms,
      num_of_iron_atoms,
      atoms,
      velocity_schedules,
    }
  }

  pub fn from_slice(atoms: &[Particle]) -> Self {
    Self::new(atoms.to_vec())
  }
}
