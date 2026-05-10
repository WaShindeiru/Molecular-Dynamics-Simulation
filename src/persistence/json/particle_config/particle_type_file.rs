use crate::particle::Particle;

#[derive(Clone, Copy, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "PascalCase")]
pub enum ParticleTypeFile {
  Atom,
  CustomPathAtom,
}

impl ParticleTypeFile {
  pub fn from_runtime(value: &Particle) -> Self {
    match value {
      Particle::Atom(_) => ParticleTypeFile::Atom,
      Particle::CustomPathAtom(_) => ParticleTypeFile::CustomPathAtom,
    }
  }
}
