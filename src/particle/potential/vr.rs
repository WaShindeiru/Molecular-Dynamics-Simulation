use nalgebra::Vector3;

pub fn vr(r_mag: f64, D0: f64, S: f64, Beta: f64, r0: f64) -> f64 {
  D0 / (S - 1.) * (-Beta * (2. * S).sqrt() * (r_mag - r0)).exp()
}

pub fn vr_gradient(r_ij_vec: &Vector3<f64>, r_ij_mag: f64, D0: f64, S: f64, Beta: f64, r0: f64) -> Vector3<f64> {
  let common = vr(r_ij_mag, D0, S, Beta, r0) * (Beta * (2. * S).sqrt());

  assert!(!r_ij_vec.x.is_nan() && !r_ij_vec.y.is_nan() && !r_ij_vec.z.is_nan());
  assert!(!common.is_nan());
  assert!(!r_ij_mag.is_nan());

  Vector3::new(
      common * r_ij_vec.x / r_ij_mag,
      common * r_ij_vec.y / r_ij_mag,
      common * r_ij_vec.z / r_ij_mag
  )
}