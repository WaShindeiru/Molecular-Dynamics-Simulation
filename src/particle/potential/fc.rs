use nalgebra::{Vector3};
use std::f64::consts::PI;
use crate::data::constants::get_constant;
use crate::data::{Constant, InteractionType};

pub fn fc(r_mag: f64, R: f64, D: f64) -> f64 {

  let R1 = (R-D).abs();
  let R2 = (R+D).abs();

  match r_mag {
    r_ if r_ <= R1 => 1.,
    r_ if r_ > R2 => 0.,
    r_ => 0.5 - 0.5 * (PI / (2. * D) * (r_ - R)).sin(),
  }
}

pub fn fc_gradient(r_ij_vec: &Vector3<f64>, r_ij_mag: f64, R: f64, D: f64) -> Vector3<f64> {
  let R1 = (R-D).abs();
  let R2 = (R+D).abs();

  match r_ij_mag {
    rij_ if rij_ <= R1 => Vector3::new(0., 0., 0.),
    rij_ if rij_ > R2 => Vector3::new(0., 0., 0.),
    rij_ => {
      let common = PI / (4.0 * D) * (PI / (2. * D) * (rij_ - R)).cos();
      Vector3::new(
        common * r_ij_vec.x / rij_,
        common * r_ij_vec.y / rij_,
        common * r_ij_vec.z / rij_)
    }
  }
}