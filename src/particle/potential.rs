use std::collections::HashMap;
use nalgebra::Vector3;
use crate::data::constants::get_constant;
use crate::data::{Constant, InteractionType};
use crate::data::types::get_interaction_type;
use crate::particle::{Atom, ParticleOperations};
use crate::particle::atom_collection::{AtomCollection, AtomMetadata};
use crate::particle::potential::b::g;
use crate::particle::potential::fc::{fc, fc_gradient};
use crate::particle::potential::va::{va, va_gradient};
use crate::particle::potential::vr::{vr, vr_gradient};
use crate::utils::math::cos_from_vec;

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

pub fn compute_forces_potential(particles: &Vec<Box<dyn ParticleOperations>>) -> FPInfo {
  let mut result: Vec<FP> = vec![
    FP {
      force: Vector3::new(0., 0., 0.),
      potential_energy: 0.,
    };
    particles.len()
  ];

  let mut gradients_cache: Vec<Vector3<f64>> = vec![Vector3::new(0., 0., 0.); particles.len()];
  let mut potential_energy_total: f64 = 0.;

  for (i, particle_i) in particles.iter().enumerate() {
    let mut neigh: Vec<usize> = Vec::new();
    for j in 0..particles.len() {
      if j != i {
        neigh.push(j);
      }
    }

    for j_ in neigh.iter() {
      let j = *j_;
      assert_ne!(j, i);

      gradients_cache = vec![Vector3::new(0., 0., 0.); particles.len()];

      let partile_j = particles.get(j).unwrap();
      let interaction_type_ij = get_interaction_type(particle_i.get_type(), partile_j.get_type());

      let mut bij_grad_i = Vector3::new(0., 0., 0.);
      let mut bij_grad_j = Vector3::new(0., 0., 0.);

      let R_ij = get_constant(&interaction_type_ij, Constant::R);
      let D_ij = get_constant(&interaction_type_ij, Constant::D);
      let S_ij = get_constant(&interaction_type_ij, Constant::S);
      let D0_ij = get_constant(&interaction_type_ij, Constant::D0);
      let Beta_ij = get_constant(&interaction_type_ij, Constant::Beta);
      let r0_ij = get_constant(&interaction_type_ij, Constant::r0);

      let r_ij_vec = partile_j.get_position() - particle_i.get_position();
      let r_ij_mag = r_ij_vec.magnitude();

      let fc_ij = fc(r_ij_mag, R_ij, D_ij);
      let va_ij = va(r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);
      let vr_ij = vr(r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);

      let fc_ij_grad_i = fc_gradient(&r_ij_vec, r_ij_mag, R_ij, D_ij);
      let vr_ij_grad_i = vr_gradient(&r_ij_vec, r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);
      let va_ij_grad_i = va_gradient(&r_ij_vec, r_ij_mag, D0_ij, S_ij, Beta_ij, r0_ij);

      let mut chi_ij: f64 = 0.;

      for _k in neigh.iter() {
        let k = *_k;
        if k == j {continue};

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

      bij_grad_i = bij_grad_i * b_ij_grad_chi_ij;
      bij_grad_j = bij_grad_j * b_ij_grad_chi_ij;

      for k_ in neigh.iter() {
        let k = *k_;
        if k == j || k == i {continue};
        gradients_cache[k] = gradients_cache[k] * b_ij_grad_chi_ij;
      }

      let force_i = -0.5 * (fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (vr_ij_grad_i - bij_grad_i * va_ij - b_ij * va_ij_grad_i));
      result[i].force += force_i;

      let force_j = -0.5 * (-fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (-vr_ij_grad_i - bij_grad_i * va_ij - b_ij * (-va_ij_grad_i)));
      result[j].force += force_j;

      for k_ in neigh.iter() {
        let k = *k_;
        if k == j || k == i {continue};
        let force_k = -0.5 * ( -fc_ij * gradients_cache[k] * va_ij );
        result[k].force += force_k;
      }

      potential_energy_total = potential_energy_total + 0.5 * fc_ij * (vr_ij - b_ij * va_ij);
    }
  }

  FPInfo {
    fp: result,
    potential_energy: potential_energy_total,
  }
}

// pub fn compute_force_i(atom_cont: &dyn AtomCollection, atom_i: &dyn AtomMetadata) -> Vector3<f64> {
//
//   let mut result = Vector3::new(0., 0., 0.);
//   let i_id = atom_i.get_id();
//
//   for (atom_j_id, atom_j) in atom_cont.get_all_atoms().iter() {
//     assert_eq!(*atom_j_id, atom_j.get_id());
//     if *atom_j_id == i_id {
//       continue;
//     }
//
//     // remove this
//     assert_eq!(*atom_i.get_position(), *atom_cont.get_atom_by_id(i_id).unwrap().get_position());
//     let r_ij_vec = atom_j.get_position() - atom_cont.get_atom_by_id(i_id).unwrap().get_position();
//     let r_ij_mag = r_ij_vec.magnitude();
//     let interaction_type_ij = get_interaction_type(atom_i.get_type(), atom_j.get_type());
//
//     let fc_ij_grad = fc::fc_gradient(&r_ij_vec, &interaction_type_ij);
//
//     let vr_ij = vr::vr(r_ij_mag, &interaction_type_ij);
//     let b_ij = b::b(atom_cont, i_id, atom_j.get_id());
//     let va_ij = va::va(r_ij_mag, &interaction_type_ij);
//
//     let fc_ij = fc::fc(r_ij_mag, &interaction_type_ij);
//
//     let vr_ij_grad = vr::vr_gradient(&r_ij_vec, &interaction_type_ij);
//     let b_ij_grad = b::b_gradient(atom_cont, i_id, atom_j.get_id());
//     let va_ij_grad = va::va_gradient(&r_ij_vec, &interaction_type_ij);
//
//     let force_ij = -0.5 * (fc_ij_grad * (vr_ij - b_ij * va_ij)
//       + fc_ij * (vr_ij_grad - b_ij * va_ij_grad - va_ij * b_ij_grad));
//     result += force_ij;
//   }
//
//   result
// }

// pub fn compute_potential_energy_i(atom_cont: &dyn AtomCollection, atom_i: &dyn AtomMetadata) -> f64 {
//   let mut result = 0.;
//   let i_id = atom_i.get_id();
//
//   for (atom_j_id, atom_j) in atom_cont.get_all_atoms().iter() {
//     assert_eq!(*atom_j_id, atom_j.get_id());
//     if *atom_j_id == i_id {
//       continue;
//     }
//
//     // remove this assert later
//     assert_eq!(*atom_i.get_position(), *atom_cont.get_atom_by_id(i_id).unwrap().get_position());
//     let r_ij_vec = atom_j.get_position() - atom_cont.get_atom_by_id(i_id).unwrap().get_position();
//     let r_ij_mag = r_ij_vec.magnitude();
//     let interaction_type_ij = get_interaction_type(atom_i.get_type(), atom_j.get_type());
//
//     let fc_ij = fc::fc(r_ij_mag, &interaction_type_ij);
//     let vr_ij = vr::vr(r_ij_mag, &interaction_type_ij);
//     let b_ij = b::b(atom_cont, i_id, atom_j.get_id());
//     let va_ij = va::va(r_ij_mag, &interaction_type_ij);
//
//     let potential_energy_ij = 0.5 * fc_ij * (vr_ij - b_ij * va_ij);
//     result += potential_energy_ij;
//   }
//
//   result
// }