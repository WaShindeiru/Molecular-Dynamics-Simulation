use crate::data::types::AtomType;
use crate::particle::Particle;
use nalgebra::Vector3;

#[derive(Clone)]
pub struct VelocityScheduleConfig {
  pub particle_velocity_manager_id: usize,
  /// Iteration-velocity pairs; must include an entry for iteration 0.
  pub velocities: Vec<(usize, Vector3<f64>)>,
}

#[derive(Clone)]
pub struct ControlVelocityScheduleConfig {
  pub control_velocity_manager_id: usize,
  /// Iteration/component_velocity/desired_velocity triples; must include an entry for iteration 0.
  pub entries: Vec<(usize, Vector3<f64>, Vector3<f64>)>,
}

#[derive(Clone)]
pub struct ParticleConfig {
  pub atoms: Vec<Particle>,
  pub num_of_atoms: usize,
  pub num_of_carbon_atoms: usize,
  pub num_of_iron_atoms: usize,
  pub velocity_schedules: Vec<VelocityScheduleConfig>,
  pub control_velocity_schedules: Vec<ControlVelocityScheduleConfig>,
}

impl ParticleConfig {
  pub fn new(atoms: Vec<Particle>) -> Self {
    Self::new_with_schedules(atoms, vec![])
  }

  pub fn new_with_schedules(
    atoms: Vec<Particle>,
    velocity_schedules: Vec<VelocityScheduleConfig>,
  ) -> Self {
    Self::new_with_all_schedules(atoms, velocity_schedules, vec![])
  }

  pub fn new_with_all_schedules(
    atoms: Vec<Particle>,
    velocity_schedules: Vec<VelocityScheduleConfig>,
    control_velocity_schedules: Vec<ControlVelocityScheduleConfig>,
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
      control_velocity_schedules,
    }
  }

  pub fn from_slice(atoms: &[Particle]) -> Self {
    Self::new(atoms.to_vec())
  }
}
