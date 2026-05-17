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

pub fn check_position_constraint_periodic_quadratic(
  velocity: Vector3<f64>,
  acceleration: Vector3<f64>,
  previous_position: Vector3<f64>,
  time_step: f64,
  mut position: Vector3<f64>,
  container_size: Vector3<f64>,
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

    position = solve_quadratic_equation_lower_boundary(
      velocity,
      acceleration,
      previous_position,
      position,
      time_step
    );

  } else if position.z > container_size.z {
    compliant = false;
    z_comp = Compliance::ExceededHigherBoundary;

    position = solve_quadratic_equation_upper_boundary(
      velocity,
      acceleration,
      previous_position,
      position,
      time_step,
      container_size,
    );
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

fn solve_quadratic_equation_lower_boundary(
  velocity: Vector3<f64>,
  acceleration: Vector3<f64>,
  previous_position: Vector3<f64>,
  new_position: Vector3<f64>,
  time_step: f64,
) -> Vector3<f64> {

  let common = ((velocity.z).powi(2) - 2.0 * acceleration.z * previous_position.z).sqrt();

  let time_solution_1 = (-velocity.z + common) / acceleration.z;
  let time_solution_2 = (-velocity.z - common) / acceleration.z;

  assert!(time_solution_1 <= 0. || time_solution_2 <= 0.);
  assert!(time_solution_1 >= 0. || time_solution_2 >= 0.);

  let correct_solution = if time_solution_1 > 0. {
    time_solution_1
  } else if time_solution_2 > 0. {
    time_solution_2
  } else {
    panic!("no correct solution!");
  };

  assert!(correct_solution <= time_step);

  let z_vel_when_hit = velocity.z + acceleration.z * correct_solution;
  let z_inverted_vel = -z_vel_when_hit;

  let remaining_time_step = time_step - correct_solution;

  let correct_z_position = acceleration.z * remaining_time_step.powi(2) / 2.0 + z_inverted_vel * remaining_time_step;

  Vector3::new(new_position.x, new_position.y, correct_z_position)
}

fn solve_quadratic_equation_upper_boundary(
  velocity: Vector3<f64>,
  acceleration: Vector3<f64>,
  previous_position: Vector3<f64>,
  new_position: Vector3<f64>,
  time_step: f64,
  world_size: Vector3<f64>,
) -> Vector3<f64> {

  let common = ((velocity.z).powi(2) - 2.0 * acceleration.z * (previous_position.z - world_size.z)).sqrt();

  let time_solution_1 = (-velocity.z + common) / acceleration.z;
  let time_solution_2 = (-velocity.z - common) / acceleration.z;

  assert!(time_solution_1 <= 0. || time_solution_2 <= 0.);
  assert!(time_solution_1 >= 0. || time_solution_2 >= 0.);

  let correct_solution = if time_solution_1 > 0. {
    time_solution_1
  } else if time_solution_2 > 0. {
    time_solution_2
  } else {
    panic!("no correct solution!");
  };

  assert!(correct_solution <= time_step);

  let z_vel_when_hit = velocity.z + acceleration.z * correct_solution;
  let z_inverted_vel = -z_vel_when_hit;

  let remaining_time_step = time_step - correct_solution;

  let correct_z_position = acceleration.z * remaining_time_step.powi(2) / 2.0 + z_inverted_vel * remaining_time_step + world_size.z;
  
  Vector3::new(new_position.x, new_position.y, correct_z_position)
}

pub fn check_position_constraint_periodic_all(
  mut position: Vector3<f64>,
  container_size: &Vector3<f64>,
) -> (Vector3<f64>, ParticleCompliance) {
  let x_comp = Compliance::Compliant;
  let y_comp = Compliance::Compliant;
  let z_comp = Compliance::Compliant;
  let compliant = true;

  position.x = (position.x + container_size.x) % container_size.x;
  position.y = (position.y + container_size.y) % container_size.y;
  position.z = (position.z + container_size.z) % container_size.z;

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
