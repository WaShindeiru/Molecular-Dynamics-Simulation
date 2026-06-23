use nalgebra::Vector3;
use std::f64::consts::PI;

#[cfg(test)]
mod tests {
  use super::*;

  /// Numerically estimates the gradient of fc w.r.t. r_i using central differences.
  ///
  /// r_ij = r_j - r_i, so perturbing r_i by +eps changes r_ij by -eps, hence the
  /// sign flip in the central difference below.
  fn numerical_fc_gradient(r_ij_vec: &Vector3<f64>, cutoff_r: f64, cutoff_d: f64) -> Vector3<f64> {
    let eps = 1e-6;
    let axes = [
      Vector3::new(eps, 0., 0.),
      Vector3::new(0., eps, 0.),
      Vector3::new(0., 0., eps),
    ];

    let components: Vec<f64> = axes
      .iter()
      .map(|e| {
        (fc((r_ij_vec - e).magnitude(), cutoff_r, cutoff_d)
          - fc((r_ij_vec + e).magnitude(), cutoff_r, cutoff_d))
          / (2. * eps)
      })
      .collect();

    Vector3::new(components[0], components[1], components[2])
  }

  // R=2.0, D=0.5 → transition region (1.5, 2.5).
  const R: f64 = 2.0;
  const D: f64 = 0.5;

  #[test]
  fn fc_gradient_matches_finite_difference_in_transition_region() {
    let tol = 1e-5;
    // Non-axis-aligned direction exercises all three gradient components at once.
    let direction = Vector3::new(3., 4., 0.).normalize();

    for r in [1.6, 1.75, 2.0, 2.25, 2.4] {
      let r_ij_vec = direction * r;
      let analytical = fc_gradient(&r_ij_vec, r_ij_vec.magnitude(), R, D);
      let numerical = numerical_fc_gradient(&r_ij_vec, R, D);
      let diff = (analytical - numerical).magnitude();
      assert!(
        diff < tol,
        "r={r}: analytical={analytical:?} numerical={numerical:?} diff={diff}"
      );
    }
  }

  #[test]
  fn fc_gradient_is_zero_outside_transition() {
    for r in [0.5, 1.0, 1.4, 2.6, 3.0, 5.0] {
      let r_ij_vec = Vector3::new(r, 0., 0.);
      let gradient = fc_gradient(&r_ij_vec, r, R, D);
      assert_eq!(
        gradient,
        Vector3::zeros(),
        "expected zero gradient outside transition at r={r}"
      );
    }
  }

  #[test]
  fn fc_gradient_vanishes_near_transition_boundaries() {
    // At R±D the analytical gradient is exactly 0 (cos(±π/2) = 0).
    // Just inside the boundary it should be very small.
    let tol = 1e-2;
    for r in [R - D + 1e-4, R + D - 1e-4] {
      let r_ij_vec = Vector3::new(r, 0., 0.);
      let gradient = fc_gradient(&r_ij_vec, r, R, D);
      assert!(
        gradient.magnitude() < tol,
        "gradient should vanish near boundary at r={r}, got magnitude={}",
        gradient.magnitude()
      );
    }
  }
}

pub fn fc(r_mag: f64, R: f64, D: f64) -> f64 {
  let R1 = (R - D).abs();
  let R2 = (R + D).abs();

  match r_mag {
    r_ if r_ <= R1 => 1.,
    r_ if r_ > R2 => 0.,
    r_ => 0.5 - 0.5 * (PI / (2. * D) * (r_ - R)).sin(),
  }
}

pub fn fc_gradient(r_ij_vec: &Vector3<f64>, r_ij_mag: f64, R: f64, D: f64) -> Vector3<f64> {
  let R1 = (R - D).abs();
  let R2 = (R + D).abs();

  match r_ij_mag {
    rij_ if rij_ <= R1 => Vector3::new(0., 0., 0.),
    rij_ if rij_ > R2 => Vector3::new(0., 0., 0.),
    rij_ => {
      let common = PI / (4.0 * D) * (PI / (2. * D) * (rij_ - R)).cos();
      Vector3::new(
        common * r_ij_vec.x / rij_,
        common * r_ij_vec.y / rij_,
        common * r_ij_vec.z / rij_,
      )
    }
  }
}
