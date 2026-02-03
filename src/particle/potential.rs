use nalgebra::Vector3;
use crate::data::constants::get_constant;
use crate::data::{Constant};
use crate::data::types::get_interaction_type;
use crate::particle::{Particle};
use crate::particle::potential::b::g;
use crate::particle::potential::fc::{fc, fc_gradient};
use crate::particle::potential::va::{va, va_gradient};
use crate::particle::potential::vr::{vr, vr_gradient};
use crate::utils::math::cos_from_vec;

const OPTIMIZATION: bool = true;

pub mod fc;
pub mod vr;
pub mod va;
pub mod b;

#[derive(Debug, PartialEq, Clone)]
pub struct FP {
  pub force: Vector3<f64>,
  pub potential_energy: f64,
}

pub struct FPInfo {
  pub fp: Vec<FP>,
  pub potential_energy: f64,
}

// TODO: Optimize cos computation
pub fn compute_forces_potential(particles: &Vec<Particle>) -> FPInfo {
  let mut result: Vec<FP> = vec![
    FP {
      force: Vector3::new(0., 0., 0.),
      potential_energy: 0.,
    };
    particles.len()
  ];

  let mut gradients_cache: Vec<Vector3<f64>> = vec![Vector3::new(0., 0., 0.); particles.len()];
  let mut potential_energy_total: f64 = 0.;
  let mut neighbours: Vec<usize> = Vec::with_capacity(particles.len() - 1);

  // let mut num_of_neighbours: i32;
  // let mut num_of_neighbours_cache: Vec<i32> = vec![0; particles.len()];

  for (i, particle_i) in particles.iter().enumerate() {
    assert_eq!(i, particle_i.get_id() as usize);
    neighbours = Vec::with_capacity(particles.len() - 1);

    if OPTIMIZATION {
      for (j, particle_j) in particles.iter().enumerate() {
        if j == i { continue }

        let r_ij_vec = particle_j.get_position() - particle_i.get_position();
        let r_ij_mag = r_ij_vec.magnitude();
        let interaction_type_ij = get_interaction_type(particle_i.get_type(), particle_j.get_type());

        let R_ij = get_constant(&interaction_type_ij, Constant::R);
        let D_ij = get_constant(&interaction_type_ij, Constant::D);
        let fc_ij = fc::fc(r_ij_mag, R_ij, D_ij);
        if fc_ij < 1e-10 { continue }
        else {
          neighbours.push(j);
        }
      }
    } else {
      for (j, _) in particles.iter().enumerate() {
        if j == i { continue }
        neighbours.push(j);
      }
    }


    for j_ in neighbours.iter() {
      let j = *j_;
      assert_ne!(j, i);

      gradients_cache = vec![Vector3::new(0., 0., 0.); particles.len()];

      let particle_j = particles.get(j).unwrap();
      let interaction_type_ij = get_interaction_type(particle_i.get_type(), particle_j.get_type());

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
      let va_ij_grad_i = va_gradient(&r_ij_vec, r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);

      let mut bij_grad_i = Vector3::new(0., 0., 0.);
      let mut bij_grad_j = Vector3::new(0., 0., 0.);

      let mut chi_ij: f64 = 0.;

      for _k in neighbours.iter() {
        let k = *_k;
        if (k == j || k == i) {continue};

        let particle_k = particles.get(k).unwrap();
        let r_ik_vec = particle_k.get_position() - particle_i.get_position();
        let r_ik_mag = r_ik_vec.magnitude();

        let interaction_type_ik = get_interaction_type(particle_i.get_type(), particle_k.get_type());
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
        let g_ik_grads = b::g_ik_gradient(&r_ij_vec, r_ij_mag, &r_ik_vec, r_ik_mag,
                                          cos_theta_ijk, gamma_ik, c_ik, d_ik, h_ik);

        bij_grad_i += fc_ik * g_ik_grads.grad_i + g_ik * fc_ik_grad_i;
        bij_grad_j += fc_ik * g_ik_grads.grad_j;
        let bij_grad_k = fc_ik * g_ik_grads.grad_k + g_ik * fc_ik_grad_k;
        gradients_cache[k] = bij_grad_k;

        chi_ij += fc_ik * g_ik;
      }

      let b_ij = (1. / (1. + chi_ij).sqrt());
      let b_ij_grad_chi_ij = -0.5 * (1. + chi_ij).powf(-1.5);
      // let b_ij = 1.;
      // let b_ij_grad_chi_ij = 0.;

      bij_grad_i = bij_grad_i * b_ij_grad_chi_ij;
      bij_grad_j = bij_grad_j * b_ij_grad_chi_ij;

      for k_ in neighbours.iter() {
        let k = *k_;
        if k == j || k == i {continue};
        gradients_cache[k] = gradients_cache[k] * b_ij_grad_chi_ij;
      }

      let force_i = -0.5 * (fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (vr_ij_grad_i - bij_grad_i * va_ij - b_ij * va_ij_grad_i));
      result[i].force += force_i;

      let force_j = -0.5 * (-fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (-vr_ij_grad_i - bij_grad_j * va_ij - b_ij * (-va_ij_grad_i)));
      result[j].force += force_j;

      for k_ in neighbours.iter() {
        let k = *k_;
        if k == j || k == i {continue};
        let force_k = -0.5 * ( -fc_ij * gradients_cache[k] * va_ij );
        result[k].force += force_k;
      }

      let potential_energy_partial = 0.5 * fc_ij * (vr_ij - b_ij * va_ij);
      potential_energy_total = potential_energy_total + potential_energy_partial;
      result[j].potential_energy += potential_energy_partial;
    }
  }

  FPInfo {
    fp: result,
    potential_energy: potential_energy_total,
  }
}