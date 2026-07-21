use crate::particle::particle::ParticleKind;
use crate::particle::Particle;

#[derive(Clone, Copy, serde::Serialize, serde::Deserialize)]
#[serde(rename_all = "PascalCase")]
pub enum ParticleTypeFile {
  Atom,
  CustomPathAtom,
  CustomVelocityAtom,
  VelocityControlledParticle,
}

impl ParticleTypeFile {
  pub fn from_runtime(value: &Particle) -> Self {
    match value {
      Particle::Atom(_) => ParticleTypeFile::Atom,
      Particle::CustomPathAtom(_) => ParticleTypeFile::CustomPathAtom,
      Particle::CustomVelocityAtom(_) => ParticleTypeFile::CustomVelocityAtom,
      Particle::VelocityControlledParticle(_) => ParticleTypeFile::VelocityControlledParticle,
    }
  }
}

impl From<ParticleKind> for ParticleTypeFile {
  fn from(kind: ParticleKind) -> Self {
    match kind {
      ParticleKind::Atom => ParticleTypeFile::Atom,
      ParticleKind::CustomPathAtom => ParticleTypeFile::CustomPathAtom,
      ParticleKind::CustomVelocityAtom => ParticleTypeFile::CustomVelocityAtom,
      ParticleKind::VelocityControlledParticle => ParticleTypeFile::VelocityControlledParticle,
    }
  }
}
