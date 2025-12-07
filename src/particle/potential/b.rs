use nalgebra::{Vector, Vector3};
use crate::data::constants::get_constant;
use crate::data::{Constant, InteractionType};
use crate::data::types::get_interaction_type;
use crate::particle::atom_collection::AtomCollection;
use crate::sim_core::simple_atom_container::SimpleAtomContainer;
use crate::particle::potential::fc::{fc, fc_gradient};
use crate::utils::math::cos_from_vec;

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

// pub fn x_value(atom_cont: &dyn AtomCollection, i_id: u64, j_id: u64) -> f64 {
//   let i_atom = atom_cont.get_atom_by_id(i_id).unwrap();
//   let j_atom = atom_cont.get_atom_by_id(j_id).unwrap();
//   let r_ij_vec = j_atom.get_position() - i_atom.get_position();
// 
//   let atoms = atom_cont.get_all_atoms();
// 
//   let mut result = 0.;
// 
//   for (k_atom_id, k_atom) in atoms.iter() {
//     assert_eq!(*k_atom_id, k_atom.get_id());
//     if *k_atom_id == i_id || *k_atom_id == j_id {
//       continue;
//     }
// 
//     let interaction_type_ik = get_interaction_type(i_atom.get_type(), k_atom.get_type());
// 
//     let r_ik_vec = k_atom.get_position() - i_atom.get_position();
//     let cos_theta_ijk = cos_from_vec(&r_ij_vec, &r_ik_vec);
// 
//     result += fc(r_ik_vec.magnitude(), &interaction_type_ik) * g(cos_theta_ijk, &interaction_type_ik);
//   }
// 
//   result
// }
// 
// fn x_gradient_one(r_ij_vec: &Vector3<f64>, r_ik_vec: &Vector3<f64>, i_ik: &InteractionType) -> Vector3<f64> {
//   let cos_theta_ijk = cos_from_vec(r_ij_vec, r_ik_vec);
// 
//   fc_gradient(r_ik_vec, i_ik) * g(cos_theta_ijk, i_ik) + fc(r_ik_vec.magnitude(), i_ik) * g_ik_gradient(r_ij_vec, r_ik_vec, cos_theta_ijk, i_ik)
// }

// pub fn x_gradient(atom_cont: &dyn AtomCollection, i_id: u64, j_id: u64) -> Vector3<f64> {
//   let i_atom = atom_cont.get_atom_by_id(i_id).unwrap();
//   let j_atom = atom_cont.get_atom_by_id(j_id).unwrap();
//   let r_ij = j_atom.get_position() - i_atom.get_position();
// 
//   let atoms = atom_cont.get_all_atoms();
// 
//   let mut result = Vector3::new(0., 0., 0.);
// 
//   for (k_atom_id, k_atom) in atoms.iter() {
//     assert_eq!(*k_atom_id, k_atom.get_id());
//     if *k_atom_id == i_id || *k_atom_id == j_id {
//       continue;
//     }
// 
//     let r_ik = k_atom.get_position() - i_atom.get_position();
//     let i_ik = get_interaction_type(i_atom.get_type(), k_atom.get_type());
// 
//     result += x_gradient_one(&r_ij, &r_ik, &i_ik);
//   }
// 
//   result
// }

// pub fn b_gradient(atom_cont: &dyn AtomCollection, i_id: u64, j_id: u64) -> Vector3<f64> {
//   let x_val_ij = x_value(atom_cont, i_id, j_id);
//   let x_grad_ij = x_gradient(atom_cont, i_id, j_id);
// 
//   -0.5 * (1. + x_val_ij).powf(-1.5) * x_grad_ij
// }
// 
// pub fn b(atom_cont: &dyn AtomCollection, i_id: u64, j_id: u64) -> f64 {
//   let x_val_ij = x_value(atom_cont, i_id, j_id);
// 
//   (1. + x_val_ij).powf(-0.5)
// }