use nalgebra::{Vector3};

pub struct GGradient {
  pub grad_i: Vector3<f64>,
  pub grad_j: Vector3<f64>,
  pub grad_k: Vector3<f64>,
}

pub fn g_ik_gradient(r_ij_vec: &Vector3<f64>, r_ij_mag: f64, r_ik_vec: &Vector3<f64>, r_ik_mag: f64,
                     cos_ijk: f64, gamma: f64, c: f64, d: f64, h: f64) -> GGradient {
  let common = (2. * gamma * c.powi(2) * (h + cos_ijk)) / (d.powi(2) + (h + cos_ijk).powi(2)).powi(2);

  let grad_j = common * (r_ik_vec / (r_ij_mag * r_ik_mag) - cos_ijk * r_ij_vec / r_ij_mag.powi(2));
  let grad_k = common * (r_ij_vec / (r_ij_mag * r_ik_mag) - cos_ijk * r_ik_vec / r_ik_mag.powi(2));
  let grad_i = -(&grad_j + &grad_k);

  GGradient {
    grad_i,
    grad_j,
    grad_k,
  }
}

pub fn g(theta_cos: f64, gamma: f64, c: f64, d: f64, h: f64) -> f64 {
  gamma * (
    1. + (c/d).powi(2) 
      - c.powi(2)
      /
      ( d.powi(2) + (h + theta_cos).powi(2) )
  )
}