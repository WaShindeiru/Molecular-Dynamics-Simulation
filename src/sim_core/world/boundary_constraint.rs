use nalgebra::Vector3;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::Compliance::Compliant;

#[derive(Copy, Clone)]
pub struct ParticleCompliance {
  pub compliant: bool,
  pub x: Compliance,
  pub y: Compliance,
  pub z: Compliance,
}

#[derive(Copy, Clone)]
pub enum Compliance {
  Compliant,
  ExceededLowerBoundary,
  ExceededHigherBoundary,
}

pub fn check_position_constraint(mut position: Vector3<f64>, container_size: &Vector3<f64>) 
  -> (Vector3<f64>, ParticleCompliance) {
  let mut x = Compliance::Compliant;
  let mut y = Compliance::Compliant;
  let mut z = Compliance::Compliant;
  let mut compliant = true;
  
  if position.x < 0.0 {
    compliant = false;
    x = Compliance::ExceededLowerBoundary;
    position.x = - position.x;
  } else if position.x > container_size.x {
    compliant = false;
    x = Compliance::ExceededHigherBoundary;
    position.x = 2. * container_size.x - position.x;
  }

  if position.y < 0.0 {
    compliant = false;
    y = Compliance::ExceededLowerBoundary;
    position.y = - position.y;
  } else if position.y > container_size.y {
    compliant = false;
    y = Compliance::ExceededHigherBoundary;
    position.y = 2. * container_size.y - position.y;
  }

  if position.z < 0. {
    compliant = false;
    z = Compliance::ExceededLowerBoundary;
    position.z = - position.z;
  } else if position.z > container_size.z {
    compliant = false;
    z = Compliance::ExceededHigherBoundary;
    position.z = 2. * container_size.z - position.z;
  }

  (position, ParticleCompliance {
    compliant,
    x,
    y,
    z,
  })
}

pub fn apply_velocity_constraint(compliance: &ParticleCompliance, mut velocity: Vector3<f64>) 
  -> Vector3<f64> {
  velocity.x = match compliance.x {
    Compliant => velocity.x,
    Compliance::ExceededLowerBoundary 
    | Compliance::ExceededHigherBoundary => - velocity.x 
  };
  
  velocity.y = match compliance.y {
    Compliant => velocity.y,
    Compliance::ExceededLowerBoundary
    | Compliance::ExceededHigherBoundary => - velocity.y
  };

  velocity.z = match compliance.z {
    Compliant => velocity.z,
    Compliance::ExceededLowerBoundary
    | Compliance::ExceededHigherBoundary => - velocity.z
  };
  
  velocity
}