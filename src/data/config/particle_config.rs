use crate::data::types::AtomType;
use crate::particle::Particle;

#[derive(Clone)]
pub struct ParticleConfig {
  pub atoms: Vec<Particle>,
  pub num_of_atoms: usize,
  pub num_of_carbon_atoms: usize,
  pub num_of_iron_atoms: usize,
}

impl ParticleConfig {
  pub fn new(atoms: Vec<Particle>) -> Self {
    let (num_of_carbon_atoms, num_of_iron_atoms) =
      atoms
        .iter()
        .fold((0usize, 0usize), |acc, atom| match atom.get_type() {
          AtomType::C => (acc.0 + 1, acc.1),
          AtomType::Fe => (acc.0, acc.1 + 1),
        });

    ParticleConfig {
      num_of_atoms: atoms.len(),
      num_of_carbon_atoms,
      num_of_iron_atoms,
      atoms,
    }
  }

  pub fn from_slice(atoms: &[Particle]) -> Self {
    Self::new(atoms.to_vec())
  }
}
