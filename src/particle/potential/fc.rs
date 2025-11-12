use nalgebra::{Vector3};
use std::f64::consts::PI;
use crate::data::constants::get_constant;
use crate::data::{Constant, InteractionType};

pub fn fc(r_mag: f64, interaction_type: &InteractionType) -> f64 {
    let D = get_constant(interaction_type, Constant::D);
    let R = get_constant(interaction_type, Constant::R);

    let R1 = (R-D).abs();
    let R2 = (R+D).abs();

    match r_mag {
        r_ if r_ <= R1 => 1.,
        r_ if r_ > R2 => 0.,
        _ => 0.5 - 0.5 * (PI / (2. * D) * (r_mag - R)).sin(),
    }
}

pub fn fc_gradient(r_ij_vec: &Vector3<f64>, interaction_type: &InteractionType) -> Vector3<f64> {
    let D = get_constant(interaction_type, Constant::D);
    let R = get_constant(interaction_type, Constant::R);

    let R1 = (R-D).abs();
    let R2 = (R+D).abs();

    let r_ij_mag = r_ij_vec.magnitude();

    match r_ij_mag {
        rij_ if rij_ <= R1 => Vector3::new(0., 0., 0.),
        rij_ if rij_ > R2 => Vector3::new(0., 0., 0.),
        _ => {
            let common = -PI / (4.0 * D) * (PI / (2. * D) * (r_ij_mag - R)).cos();
            Vector3::new(
                - common * r_ij_vec.x / r_ij_mag,
                - common * r_ij_vec.y / r_ij_mag,
                - common * r_ij_vec.z / r_ij_mag)
        }
    }
}