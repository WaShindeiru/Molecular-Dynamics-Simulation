use crate::sim_core::world::boundary_constraint::Compliance::Compliant;
use crate::sim_core::world::boundary_constraint::{Compliance, ParticleCompliance};
use nalgebra::Vector3;

pub fn check_position_constraint_periodic(
  mut position: Vector3<f64>,
  container_size: &Vector3<f64>,
) -> (Vector3<f64>, ParticleCompliance) {
  let mut x_comp = Compliance::Compliant;
  let mut y_comp = Compliance::Compliant;
  let mut z_comp = Compliance::Compliant;
  let mut compliant = true;

  position.x = (position.x + container_size.x) % container_size.x;
  position.y = (position.y + container_size.y) % container_size.y;

  if position.z < 0. {
    compliant = false;
    z_comp = Compliance::ExceededLowerBoundary;
    position.z = -position.z;
  } else if position.z > container_size.z {
    compliant = false;
    z_comp = Compliance::ExceededHigherBoundary;
    position.z = 2. * container_size.z - position.z;
  }

  (
    position,
    ParticleCompliance {
      compliant,
      x: x_comp,
      y: y_comp,
      z: z_comp,
    },
  )
}

pub fn apply_velocity_constraint_periodic(
  compliance: &ParticleCompliance,
  mut velocity: Vector3<f64>,
) -> Vector3<f64> {
  velocity.z = match compliance.z {
    Compliant => velocity.z,
    Compliance::ExceededLowerBoundary | Compliance::ExceededHigherBoundary => -velocity.z,
  };

  velocity
}

#[cfg(test)]
mod tests {
  use super::*;

  #[test]
  fn check_position_constraint_periodic_keeps_position_inside_bounds() {
    let container_size = Vector3::new(10.0, 20.0, 30.0);
    let position = Vector3::new(2.0, 3.0, 4.0);

    let (validated_position, compliance) =
      check_position_constraint_periodic(position, &container_size);

    assert_eq!(validated_position, Vector3::new(2.0, 3.0, 4.0));
    assert_eq!(
      compliance,
      ParticleCompliance {
        compliant: true,
        x: Compliance::Compliant,
        y: Compliance::Compliant,
        z: Compliance::Compliant,
      }
    );
  }

  #[test]
  fn check_position_constraint_periodic_wraps_x_and_y() {
    let container_size = Vector3::new(10.0, 20.0, 30.0);
    let position = Vector3::new(12.0, -3.0, 4.0);

    let (validated_position, compliance) =
      check_position_constraint_periodic(position, &container_size);

    assert_eq!(validated_position, Vector3::new(2.0, 17.0, 4.0));
    assert_eq!(
      compliance,
      ParticleCompliance {
        compliant: true,
        x: Compliance::Compliant,
        y: Compliance::Compliant,
        z: Compliance::Compliant,
      }
    );
  }

  #[test]
  fn check_position_constraint_periodic_wraps_x_and_x() {
    let container_size = Vector3::new(10.0, 20.0, 30.0);
    let position = Vector3::new(-2.0, 5.0, 4.0);

    let (validated_position, compliance) =
      check_position_constraint_periodic(position, &container_size);

    assert_eq!(validated_position, Vector3::new(8.0, 5.0, 4.0));
    assert_eq!(
      compliance,
      ParticleCompliance {
        compliant: true,
        x: Compliance::Compliant,
        y: Compliance::Compliant,
        z: Compliance::Compliant,
      }
    );
  }

  #[test]
  fn check_position_constraint_periodic_reflects_lower_z_boundary() {
    let container_size = Vector3::new(10.0, 20.0, 30.0);
    let position = Vector3::new(2.0, 3.0, -4.0);

    let (validated_position, compliance) =
      check_position_constraint_periodic(position, &container_size);

    assert_eq!(validated_position, Vector3::new(2.0, 3.0, 4.0));
    assert_eq!(
      compliance,
      ParticleCompliance {
        compliant: false,
        x: Compliance::Compliant,
        y: Compliance::Compliant,
        z: Compliance::ExceededLowerBoundary,
      }
    );
  }

  #[test]
  fn check_position_constraint_periodic_reflects_upper_z_boundary() {
    let container_size = Vector3::new(10.0, 20.0, 30.0);
    let position = Vector3::new(2.0, 3.0, 34.0);

    let (validated_position, compliance) =
      check_position_constraint_periodic(position, &container_size);

    assert_eq!(validated_position, Vector3::new(2.0, 3.0, 26.0));
    assert_eq!(
      compliance,
      ParticleCompliance {
        compliant: false,
        x: Compliance::Compliant,
        y: Compliance::Compliant,
        z: Compliance::ExceededHigherBoundary,
      }
    );
  }
}
