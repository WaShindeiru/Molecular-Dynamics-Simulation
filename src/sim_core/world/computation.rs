use crate::data::constants::get_constants;
use crate::data::types::{AtomType, get_interaction_type};
use crate::particle::potential::b::g;
use crate::particle::potential::fc::{fc, fc_gradient};
use crate::particle::potential::va::{va, va_gradient};
use crate::particle::potential::vr::{vr, vr_gradient};
use crate::particle::potential::{b, fc};
use crate::utils::math::cos_from_vec;
use nalgebra::Vector3;
use std::any::Any;

#[derive(Debug, PartialEq, Clone)]
pub struct FP {
  pub force: Vector3<f64>,
  pub potential_energy: f64,
}

pub fn compute_new_velocity(
  half_velocity: Vector3<f64>,
  acceleration: Vector3<f64>,
  thermostat_epsilon: f64,
  time_step: f64,
) -> Vector3<f64> {
  let numerator = half_velocity + 0.5 * acceleration * time_step;
  let denominator = 1.0 + 0.5 * time_step * thermostat_epsilon;
  numerator / denominator
}

pub trait ForceComputationOperations {
  fn as_any(&self) -> &dyn Any;
  fn get_id(&self) -> usize;
  fn get_position(&self) -> Vector3<f64>;
  fn get_type(&self) -> AtomType;
  fn get_mass(&self) -> f64;
  fn prototype_clone(&self) -> Box<dyn ForceComputationOperations>;
}

impl ForceComputationOperations for Box<dyn ForceComputationOperations> {
  fn as_any(&self) -> &dyn Any { self.as_ref().as_any() }
  fn get_id(&self) -> usize { self.as_ref().get_id() }
  fn get_position(&self) -> Vector3<f64> { self.as_ref().get_position() }
  fn get_type(&self) -> AtomType { self.as_ref().get_type() }
  fn get_mass(&self) -> f64 { self.as_ref().get_mass() }
  fn prototype_clone(&self) -> Box<dyn ForceComputationOperations> { self.as_ref().prototype_clone() }
}

impl<T: ForceComputationOperations + ?Sized> ForceComputationOperations for &T {
  fn as_any(&self) -> &dyn Any { (*self).as_any() }
  fn get_id(&self) -> usize { (*self).get_id() }
  fn get_position(&self) -> Vector3<f64> { (*self).get_position() }
  fn get_type(&self) -> AtomType { (*self).get_type() }
  fn get_mass(&self) -> f64 { (*self).get_mass() }
  fn prototype_clone(&self) -> Box<dyn ForceComputationOperations> { (*self).prototype_clone() }
}

pub fn compute_forces_potential<I>(
  particles_i: I,
  particles_j: I,
  fp: &mut Vec<FP>,
  gradients_cache: &mut Vec<Vector3<f64>>,
) -> f64
where
  I: IntoIterator + Clone,
  I::Item: ForceComputationOperations,
{
  let zero_fp = FP { force: Vector3::zeros(), potential_energy: 0. };

  for temp_j in particles_j.clone().into_iter() {
    fp[temp_j.get_id()] = zero_fp.clone();
  }

  let mut potential_energy_total: f64 = 0.;

  for temp_i in particles_i.into_iter() {
    let i_id = temp_i.get_id();

    for temp_j in particles_j.clone().into_iter() {
      let j_id = temp_j.get_id();
      if j_id == i_id {
        continue;
      }

      let c_ij = get_constants(&get_interaction_type(&temp_i.get_type(), &temp_j.get_type()));

      let r_ij_vec = temp_j.get_position() - temp_i.get_position();
      let r_ij_mag = r_ij_vec.magnitude();

      let fc_ij = fc(r_ij_mag, c_ij.R, c_ij.D);
      let va_ij = va(r_ij_mag, c_ij.D0, c_ij.S, c_ij.Beta, c_ij.r0);
      let vr_ij = vr(r_ij_mag, c_ij.D0, c_ij.S, c_ij.Beta, c_ij.r0);

      let fc_ij_grad_i = fc_gradient(&r_ij_vec, r_ij_mag, c_ij.R, c_ij.D);
      let vr_ij_grad_i = vr_gradient(&r_ij_vec, r_ij_mag, c_ij.D0, c_ij.S, c_ij.Beta, c_ij.r0);
      assert!(!vr_ij_grad_i.x.is_nan() && !vr_ij_grad_i.y.is_nan() && !vr_ij_grad_i.z.is_nan());
      let va_ij_grad_i = va_gradient(&r_ij_vec, r_ij_mag, c_ij.D0, c_ij.S, c_ij.Beta, c_ij.r0);

      let mut bij_grad_i: Vector3<f64> = Vector3::zeros();
      let mut bij_grad_j: Vector3<f64> = Vector3::zeros();
      let mut chi_ij: f64 = 0.;

      for temp_k in particles_j.clone().into_iter() {
        let k_id = temp_k.get_id();
        if k_id == j_id || k_id == i_id {
          continue;
        }

        let c_ik = get_constants(&get_interaction_type(&temp_i.get_type(), &temp_k.get_type()));

        let r_ik_vec = temp_k.get_position() - temp_i.get_position();
        let r_ik_mag = r_ik_vec.magnitude();

        let fc_ik = fc(r_ik_mag, c_ik.R, c_ik.D);
        let cos_theta_ijk = cos_from_vec(&r_ij_vec, &r_ik_vec);
        let g_ik = g(cos_theta_ijk, c_ik.Gamma, c_ik.c, c_ik.d, c_ik.h);
        let fc_ik_grad_i = fc_gradient(&r_ik_vec, r_ik_mag, c_ik.R, c_ik.D);
        let fc_ik_grad_k = -fc_ik_grad_i;
        let g_ik_grads = b::g_ik_gradient(
          &r_ij_vec,
          r_ij_mag,
          &r_ik_vec,
          r_ik_mag,
          cos_theta_ijk,
          c_ik.Gamma,
          c_ik.c,
          c_ik.d,
          c_ik.h,
        );

        bij_grad_i += fc_ik * g_ik_grads.grad_i + g_ik * fc_ik_grad_i;
        bij_grad_j += fc_ik * g_ik_grads.grad_j;
        gradients_cache[k_id] = fc_ik * g_ik_grads.grad_k + g_ik * fc_ik_grad_k;

        chi_ij += fc_ik * g_ik;
      }

      let b_ij = 1. / (1. + chi_ij).sqrt();
      let b_ij_grad_chi_ij = -0.5 * (1. + chi_ij).powf(-1.5);

      bij_grad_i *= b_ij_grad_chi_ij;
      bij_grad_j *= b_ij_grad_chi_ij;

      let force_i = (fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (vr_ij_grad_i - bij_grad_i * va_ij - b_ij * va_ij_grad_i))
        * -0.5;
      assert!(!force_i.x.is_nan() && !force_i.y.is_nan() && !force_i.z.is_nan());
      fp[i_id].force += force_i;

      let force_j = (-fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (-vr_ij_grad_i - bij_grad_j * va_ij - b_ij * (-va_ij_grad_i)))
        * -0.5;
      assert!(!force_j.x.is_nan() && !force_j.y.is_nan() && !force_j.z.is_nan());
      fp[j_id].force += force_j;

      for temp_k in particles_j.clone().into_iter() {
        let k_id = temp_k.get_id();
        if k_id == j_id || k_id == i_id {
          continue;
        }
        let force_k = (-fc_ij * gradients_cache[k_id] * b_ij_grad_chi_ij * va_ij) * -0.5;
        assert!(!force_k.x.is_nan() && !force_k.y.is_nan() && !force_k.z.is_nan());
        fp[k_id].force += force_k;
      }

      let potential_energy_partial = 0.5 * fc_ij * (vr_ij - b_ij * va_ij);
      potential_energy_total += potential_energy_partial;
      fp[i_id].potential_energy += potential_energy_partial;
    }
  }

  potential_energy_total
}
