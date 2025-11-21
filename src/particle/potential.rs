use nalgebra::Vector3;
use crate::data::types::get_interaction_type;
use crate::particle::Atom;
use crate::particle::atom_collection::{AtomCollection, AtomMetadata};

pub mod fc;
pub mod vr;
pub mod va;
pub mod b;

pub fn compute_force_i(atom_cont: &dyn AtomCollection, atom_i: &dyn AtomMetadata) -> Vector3<f64> {

  let mut result = Vector3::new(0., 0., 0.);
  let i_id = atom_i.get_id();

  for (atom_j_id, atom_j) in atom_cont.get_all_atoms().iter() {
    assert_eq!(*atom_j_id, atom_j.get_id());
    if *atom_j_id == i_id {
      continue;
    }

    // remove this
    assert_eq!(*atom_i.get_position(), *atom_cont.get_atom_by_id(i_id).unwrap().get_position());
    let r_ij_vec = atom_j.get_position() - atom_cont.get_atom_by_id(i_id).unwrap().get_position();
    let r_ij_mag = r_ij_vec.magnitude();
    let interaction_type_ij = get_interaction_type(atom_i.get_type(), atom_j.get_type());

    let fc_ij_grad = fc::fc_gradient(&r_ij_vec, &interaction_type_ij);

    let vr_ij = vr::vr(r_ij_mag, &interaction_type_ij);
    let b_ij = b::b(atom_cont, i_id, atom_j.get_id());
    let va_ij = va::va(r_ij_mag, &interaction_type_ij);

    let fc_ij = fc::fc(r_ij_mag, &interaction_type_ij);

    let vr_ij_grad = vr::vr_gradient(&r_ij_vec, &interaction_type_ij);
    let b_ij_grad = b::b_gradient(atom_cont, i_id, atom_j.get_id());
    let va_ij_grad = va::va_gradient(&r_ij_vec, &interaction_type_ij);

    let force_ij = -0.5 * (fc_ij_grad * (vr_ij - b_ij * va_ij)
      + fc_ij * (vr_ij_grad - b_ij * va_ij_grad - va_ij * b_ij_grad));
    result += force_ij;
  }

  result
}

pub fn compute_potential_energy_i(atom_cont: &dyn AtomCollection, atom_i: &dyn AtomMetadata) -> f64 {
  let mut result = 0.;
  let i_id = atom_i.get_id();
  
  for (atom_j_id, atom_j) in atom_cont.get_all_atoms().iter() {
    assert_eq!(*atom_j_id, atom_j.get_id());
    if *atom_j_id == i_id {
      continue;
    }

    // remove this assert later
    assert_eq!(*atom_i.get_position(), *atom_cont.get_atom_by_id(i_id).unwrap().get_position());
    let r_ij_vec = atom_j.get_position() - atom_cont.get_atom_by_id(i_id).unwrap().get_position();
    let r_ij_mag = r_ij_vec.magnitude();
    let interaction_type_ij = get_interaction_type(atom_i.get_type(), atom_j.get_type());
    
    let fc_ij = fc::fc(r_ij_mag, &interaction_type_ij);
    let vr_ij = vr::vr(r_ij_mag, &interaction_type_ij);
    let b_ij = b::b(atom_cont, i_id, atom_j.get_id());
    let va_ij = va::va(r_ij_mag, &interaction_type_ij);
    
    let potential_energy_ij = 0.5 * fc_ij * (vr_ij - b_ij * va_ij);
    result += potential_energy_ij;
  }
  
  result
}