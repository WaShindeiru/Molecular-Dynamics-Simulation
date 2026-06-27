use crate::data::Constant;
use crate::data::constants::get_constant;
use crate::data::types::{AtomType, get_interaction_type};
use crate::particle::potential::b::g;
use crate::particle::potential::fc::{fc, fc_gradient};
use crate::particle::potential::va::{va, va_gradient};
use crate::particle::potential::vr::{vr, vr_gradient};
use crate::particle::potential::{b, fc};
use crate::utils::math::cos_from_vec;
use nalgebra::Vector3;
use std::any::Any;
use std::collections::HashMap;

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

#[derive(Debug, PartialEq, Clone)]
pub struct FP {
  pub force: Vector3<f64>,
  pub potential_energy: f64,
}

pub struct FPInfoBoxed {
  pub fp: HashMap<usize, FP>,
  pub potential_energy: f64,
  pub optimization_considered: usize,
  pub optimization_ignored: usize,
}

pub trait ForceComputationOperations {
  fn as_any(&self) -> &dyn Any;
  fn get_id(&self) -> usize;
  fn get_position(&self) -> Vector3<f64>;
  fn get_type(&self) -> AtomType;
  fn get_mass(&self) -> f64;
  fn prototype_clone(&self) -> Box<dyn ForceComputationOperations>;
}

pub fn compute_forces_potential<I>(particles_i: I, particles_j: I, optimization: bool) -> FPInfoBoxed
where
  I: IntoIterator + Clone,
  I::Item: AsRef<dyn ForceComputationOperations>,
{
  let mut fp: HashMap<usize, FP> = HashMap::new();

  let defaultFP = || FP {
    force: Vector3::new(0., 0., 0.),
    potential_energy: 0.,
  };

  for temp_j in particles_j.clone().into_iter() {
    let particle_j = temp_j.as_ref();
    let j_id = particle_j.get_id();
    fp.insert(j_id, defaultFP());
  }

  let mut optimization_considered: usize = 0;
  let mut optimization_ignored: usize = 0;

  let mut gradients_cache: HashMap<usize, Vector3<f64>> = HashMap::new();
  let mut potential_energy_total: f64 = 0.;
  let mut neighbours: Vec<usize> = Vec::new();

  let mut particles_j_cache: HashMap<usize, Box<dyn ForceComputationOperations>> = HashMap::new();

  for temp_i in particles_i.into_iter() {
    let particle_i = temp_i.as_ref();
    let i_id = particle_i.get_id();
    neighbours.clear();

    if optimization {
      for temp_j in particles_j.clone().into_iter() {
        let particle_j = temp_j.as_ref();
        let j_id = particle_j.get_id();
        if i_id == j_id {
          continue;
        }

        let r_ij_vec = particle_j.get_position() - particle_i.get_position();
        let r_ij_mag = r_ij_vec.magnitude();
        let interaction_type_ij =
          get_interaction_type(&particle_i.get_type(), &particle_j.get_type());

        let R_ij = get_constant(&interaction_type_ij, Constant::R);
        let D_ij = get_constant(&interaction_type_ij, Constant::D);
        let fc_ij = fc::fc(r_ij_mag, R_ij, D_ij);
        if fc_ij < 1e-10 {
          optimization_ignored += 1;
          continue;
        } else {
          optimization_considered += 1;
          neighbours.push(j_id);
          particles_j_cache.insert(j_id, particle_j.prototype_clone());
        }
      }
    } else {
      for temp_j in particles_j.clone().into_iter() {
        optimization_considered += 1;
        let particle_j = temp_j.as_ref();
        let j_id = particle_j.get_id();
        if i_id == j_id {
          continue;
        }
        neighbours.push(j_id);
        particles_j_cache.insert(j_id, particle_j.prototype_clone());
      }
    }

    for j_id_ in neighbours.iter() {
      let j_id = *j_id_;
      assert_ne!(i_id, j_id);

      gradients_cache.clear();

      let particle_j = particles_j_cache.get(j_id_).unwrap();
      let interaction_type_ij =
        get_interaction_type(&particle_i.get_type(), &particle_j.get_type());

      let R_ij = get_constant(&interaction_type_ij, Constant::R);
      let D_ij = get_constant(&interaction_type_ij, Constant::D);
      let S_ij = get_constant(&interaction_type_ij, Constant::S);
      let D0_ij = get_constant(&interaction_type_ij, Constant::D0);
      let Beta_ij = get_constant(&interaction_type_ij, Constant::Beta);
      let r0_ij = get_constant(&interaction_type_ij, Constant::r0);

      let r_ij_vec = particle_j.get_position() - particle_i.get_position();
      let r_ij_mag = r_ij_vec.magnitude();

      let fc_ij = fc(r_ij_mag, R_ij, D_ij);
      let va_ij = va(r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);
      let vr_ij = vr(r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);

      let fc_ij_grad_i = fc_gradient(&r_ij_vec, r_ij_mag, R_ij, D_ij);
      let vr_ij_grad_i = vr_gradient(&r_ij_vec, r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);
      assert!(!vr_ij_grad_i.x.is_nan() && !vr_ij_grad_i.y.is_nan() && !vr_ij_grad_i.z.is_nan());
      let va_ij_grad_i = va_gradient(&r_ij_vec, r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);

      let mut bij_grad_i: Vector3<f64> = Vector3::new(0., 0., 0.);
      let mut bij_grad_j: Vector3<f64> = Vector3::new(0., 0., 0.);

      let mut chi_ij: f64 = 0.;

      for k_id_ in neighbours.iter() {
        let k_id = *k_id_;
        if k_id == j_id || k_id == i_id {
          continue;
        }

        let particle_k = particles_j_cache.get(k_id_).unwrap();
        let r_ik_vec = particle_k.get_position() - particle_i.get_position();
        let r_ik_mag = r_ik_vec.magnitude();

        let interaction_type_ik =
          get_interaction_type(&particle_i.get_type(), &particle_k.get_type());
        let R_ik = get_constant(&interaction_type_ik, Constant::R);
        let D_ik = get_constant(&interaction_type_ik, Constant::D);
        let gamma_ik = get_constant(&interaction_type_ik, Constant::Gamma);
        let c_ik = get_constant(&interaction_type_ik, Constant::c);
        let d_ik = get_constant(&interaction_type_ik, Constant::d);
        let h_ik = get_constant(&interaction_type_ik, Constant::h);

        let fc_ik = fc(r_ik_mag, R_ik, D_ik);
        let cos_theta_ijk = cos_from_vec(&r_ij_vec, &r_ik_vec);
        let g_ik = g(cos_theta_ijk, gamma_ik, c_ik, d_ik, h_ik);
        let fc_ik_grad_i = fc_gradient(&r_ik_vec, r_ik_mag, R_ik, D_ik);
        let fc_ik_grad_k = -&fc_ik_grad_i;
        let g_ik_grads = b::g_ik_gradient(
          &r_ij_vec,
          r_ij_mag,
          &r_ik_vec,
          r_ik_mag,
          cos_theta_ijk,
          gamma_ik,
          c_ik,
          d_ik,
          h_ik,
        );

        bij_grad_i += fc_ik * g_ik_grads.grad_i + g_ik * fc_ik_grad_i;
        bij_grad_j += fc_ik * g_ik_grads.grad_j;
        let bij_grad_k = fc_ik * g_ik_grads.grad_k + g_ik * fc_ik_grad_k;
        gradients_cache.insert(k_id, bij_grad_k);

        chi_ij += fc_ik * g_ik;
      }

      let b_ij = 1. / (1. + chi_ij).sqrt();
      let b_ij_grad_chi_ij = -0.5 * (1. + chi_ij).powf(-1.5);

      bij_grad_i = bij_grad_i * b_ij_grad_chi_ij;
      bij_grad_j = bij_grad_j * b_ij_grad_chi_ij;

      for k_id_ in neighbours.iter() {
        let k_id = *k_id_;
        if k_id == j_id || k_id == i_id {
          continue;
        }
        *gradients_cache.get_mut(k_id_).unwrap() =
          *gradients_cache.get_mut(k_id_).unwrap() * b_ij_grad_chi_ij;
      }

      let force_i: Vector3<f64> = (fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (vr_ij_grad_i - bij_grad_i * va_ij - b_ij * va_ij_grad_i))
        * -0.5;
      let fp_i = fp.get_mut(&i_id).unwrap();
      assert!(!force_i.x.is_nan() && !force_i.y.is_nan() && !force_i.z.is_nan());
      fp_i.force += force_i;

      let force_j = (-fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (-vr_ij_grad_i - bij_grad_j * va_ij - b_ij * (-va_ij_grad_i)))
        * -0.5;
      let fp_j = fp.get_mut(&j_id).unwrap();
      assert!(!force_j.x.is_nan() && !force_j.y.is_nan() && !force_j.z.is_nan());
      fp_j.force += force_j;

      for k_id_ in neighbours.iter() {
        let k_id = *k_id_;
        if k_id == j_id || k_id == i_id {
          continue;
        };
        let force_k = (-fc_ij * gradients_cache.get(k_id_).unwrap() * va_ij) * -0.5;
        let fp_k = fp.get_mut(&k_id).unwrap();
        assert!(!force_k.x.is_nan() && !force_k.y.is_nan() && !force_k.z.is_nan());
        fp_k.force += force_k;
      }

      let potential_energy_partial = 0.5 * fc_ij * (vr_ij - b_ij * va_ij);
      potential_energy_total = potential_energy_total + potential_energy_partial;
      let fp_i = fp.get_mut(&i_id).unwrap();
      fp_i.potential_energy += potential_energy_partial;
    }
  }

  FPInfoBoxed {
    fp,
    potential_energy: potential_energy_total,
    optimization_considered,
    optimization_ignored,
  }
}
