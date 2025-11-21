use nalgebra::Vector3;
use crate::data::constants::get_constant;
use crate::data::{Constant, InteractionType};

pub fn va(r_mag: f64, interaction_type: &InteractionType) -> f64 {
  let S = get_constant(interaction_type, Constant::S);
  let D0 = get_constant(interaction_type, Constant::D0);
  let Beta = get_constant(interaction_type, Constant::Beta);
  let r0 = get_constant(interaction_type, Constant::r0);

  (S * D0) / (S - 1.) * (-Beta * (2. / S).sqrt() * (r_mag - r0)).exp()
}

pub fn va_gradient(r_ij_vec: &Vector3<f64>, interaction_type: &InteractionType) -> Vector3<f64> {
  let Beta = get_constant(interaction_type, Constant::Beta);
  let S = get_constant(interaction_type, Constant::S);

  let r_ij_mag = r_ij_vec.magnitude();
  let common = va(r_ij_mag, interaction_type) * (Beta * (2. / S).sqrt());

  Vector3::new(
    common * r_ij_vec.x / r_ij_mag,
    common * r_ij_vec.y / r_ij_mag,
    common * r_ij_vec.z / r_ij_mag
  )
}