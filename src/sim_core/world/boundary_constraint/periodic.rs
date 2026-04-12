use nalgebra::Vector3;
use crate::sim_core::world::boundary_constraint::{Compliance, ParticleCompliance};
use crate::sim_core::world::boundary_constraint::Compliance::Compliant;

pub fn check_position_constraint_periodic(mut position: Vector3<f64>, container_size: &Vector3<f64>)
                                          -> (Vector3<f64>, ParticleCompliance) {
  let mut x_comp = Compliance::Compliant;
  let mut y_comp = Compliance::Compliant;
  let mut z_comp = Compliance::Compliant;
  let mut compliant = true;

  position.x = (position.x + container_size.x) % container_size.x;
  position.y = (position.y + container_size.y) % container_size.y;

  if position.z < 0. {
    compliant = false;
    z_comp = Compliance::ExceededLowerBoundary;
    position.z = - position.z;
  } else if position.z > container_size.z {
    compliant = false;
    z_comp = Compliance::ExceededHigherBoundary;
    position.z = 2. * container_size.z - position.z;
  }

  (position, ParticleCompliance {
    compliant,
    x: x_comp,
    y: y_comp,
    z: z_comp,
  })
}

pub fn apply_velocity_constraint_periodic(compliance: &ParticleCompliance, mut velocity: Vector3<f64>) 
  -> Vector3<f64> {
  
  velocity.z = match compliance.z {
    Compliant => velocity.z,
    Compliance::ExceededLowerBoundary
    | Compliance::ExceededHigherBoundary => - velocity.z
  };

  velocity
}