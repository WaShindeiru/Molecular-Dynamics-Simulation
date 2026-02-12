use std::collections::HashMap;
use nalgebra::Vector3;
use crate::data::Constant;
use crate::data::constants::get_constant;
use crate::data::types::get_interaction_type;
use crate::particle::potential::b::g;
use crate::particle::potential::{b, fc};
use crate::particle::potential::fc::{fc, fc_gradient};
use crate::particle::potential::va::{va, va_gradient};
use crate::particle::potential::vr::{vr, vr_gradient};
use crate::particle::SimpleAtomContainer;
use crate::utils::math::cos_from_vec;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_task::BoxTask;
use crate::sim_core::world::boxed_world::BoxedWorld;

const OPTIMIZATION: bool = true;

pub fn verlet_noose_hoover_half_velocity_position<I>(previous_atom_container: I, time_step: f64,
                                                     previous_thermostat_epsilon: f64,
                                                     atom_count: usize, current_iteration: usize)
  -> (HashMap<usize, Vector3<f64>>, HashMap<usize, Particle>)
where
  I: IntoIterator,
  I::Item: AsRef<Particle>,
{
  let mut half_velocity_cache: HashMap<usize, Vector3<f64>> = HashMap::with_capacity(atom_count);
  let mut new_position_atoms: HashMap<usize, Particle> = HashMap::with_capacity(atom_count);

  for temp_i in previous_atom_container.into_iter() {
    let atom_i = temp_i.as_ref();
    let i_id = atom_i.get_id() as usize;

    let thermostat_difference = atom_i.get_acceleration() -
      previous_thermostat_epsilon * atom_i.get_velocity();
    let half_velocity_i: Vector3<f64> = atom_i.get_velocity() + thermostat_difference *
      (time_step / 2.0);
    half_velocity_cache.insert(i_id, half_velocity_i);

    let previous_position = atom_i.get_position();
    let next_position: Vector3<f64> = previous_position + half_velocity_i * time_step;
    let thermostat_work;

    if current_iteration == 0 {
      thermostat_work = 0.;
    } else {
      let thermostat_force = previous_thermostat_epsilon * atom_i.get_mass() *
        atom_i.get_velocity();
      let thermostat_path = next_position - previous_position;
      thermostat_work = thermostat_force.magnitude() * thermostat_path.magnitude()
        * cos_from_vec(&thermostat_force, &thermostat_path);
    }

    let mut new_atom_data = atom_i.custom_clone();
    new_atom_data.update_position(next_position);
    new_atom_data.set_thermostat_work(thermostat_work);

    new_position_atoms.insert(i_id, new_atom_data);
  }

  (half_velocity_cache, new_position_atoms)
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

pub fn compute_forces_potential<I>(particles_i: &I, particles_j: &I, particles_j_count: usize)
  -> FPInfoBoxed
where
  I: ?Sized,
  for<'a> &'a I: IntoIterator,
  for<'a> <&'a I as IntoIterator>::Item: AsRef<Particle>,
{
  let mut fp: HashMap<usize, FP> = HashMap::new();
  let mut optimization_considered: usize = 0;
  let mut optimization_ignored: usize = 0;

  let mut gradients_cache: HashMap<usize, Vector3<f64>> = HashMap::new();
  let mut potential_energy_total: f64 = 0.;
  let mut neighbours: Vec<usize> = Vec::with_capacity(particles_j_count);

  let mut particles_j_cache: HashMap<usize, Particle> = HashMap::new();

  for temp_i in particles_i.into_iter() {
    let particle_i = temp_i.as_ref();
    let i_id = particle_i.get_id() as usize;
    neighbours = Vec::with_capacity(particles_j_count);

    if OPTIMIZATION {
      for temp_j in particles_j.into_iter() {
        let particle_j = temp_j.as_ref();
        let j_id = particle_j.get_id() as usize;
        if i_id == j_id {
          continue
        }

        let r_ij_vec = particle_j.get_position() - particle_i.get_position();
        let r_ij_mag = r_ij_vec.magnitude();
        let interaction_type_ij = get_interaction_type(&particle_i.get_type(),
                                                       &particle_j.get_type());

        let R_ij = get_constant(&interaction_type_ij, Constant::R);
        let D_ij = get_constant(&interaction_type_ij, Constant::D);
        let fc_ij = fc::fc(r_ij_mag, R_ij, D_ij);
        if fc_ij < 1e-10 {
          optimization_ignored += 1;
          continue
        }
        else {
          optimization_considered += 1;
          neighbours.push(j_id);
          particles_j_cache.insert(j_id, particle_j.clone());
        }
      }
    } else {
      for temp_j in particles_j.into_iter() {
        optimization_considered += 1;
        let particle_j = temp_j.as_ref();
        let j_id = particle_j.get_id() as usize;
        if i_id == j_id { continue }
        neighbours.push(j_id);
        particles_j_cache.insert(j_id, particle_j.clone());
      }
    }
    
    for j_id_ in neighbours.iter() {
      let j_id = *j_id_;
      assert_ne!(i_id, j_id);

      gradients_cache = HashMap::new();

      let particle_j = particles_j_cache.get(j_id_).unwrap();
      let interaction_type_ij = get_interaction_type(&particle_i.get_type(),
                                                     &particle_j.get_type());

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

      let mut bij_grad_i: Vector3<f64> = Vector3::new(0., 0., 0.);
      let mut bij_grad_j: Vector3<f64> = Vector3::new(0., 0., 0.);

      let mut chi_ij: f64 = 0.;

      for k_id_ in neighbours.iter() {
        let k_id = *k_id_;
        if (k_id == j_id || k_id == i_id) {
          continue;
        }

        let particle_k = particles_j_cache.get(k_id_).unwrap();
        let r_ik_vec = particle_k.get_position() - particle_i.get_position();
        let r_ik_mag = r_ik_vec.magnitude();

        let interaction_type_ik = get_interaction_type(&particle_i.get_type(),
                                                       &particle_k.get_type());
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

      let defaultFP = || {
       FP {
         force: Vector3::new(0., 0., 0.),
         potential_energy: 0.
       }
      };

      let force_i: Vector3<f64> = (fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (vr_ij_grad_i - bij_grad_i * va_ij - b_ij * va_ij_grad_i)) * -0.5;
      let fp_i = fp.entry(i_id).or_insert_with(defaultFP);
      fp_i.force += force_i;

      let force_j = (-fc_ij_grad_i * (vr_ij - b_ij * va_ij)
        + fc_ij * (-vr_ij_grad_i - bij_grad_j * va_ij - b_ij * (-va_ij_grad_i))) * -0.5;
      let fp_j = fp.entry(j_id).or_insert_with(defaultFP);
      fp_j.force += force_j;

      for k_id_ in neighbours.iter() {
        let k_id = *k_id_;
        if k_id == j_id || k_id == i_id {continue};
        let force_k = (-fc_ij * gradients_cache.get(k_id_).unwrap() * va_ij) * -0.5;
        let fp_k = fp.entry(k_id).or_insert_with(defaultFP);
        fp_k.force += force_k;
      }

      let potential_energy_partial = 0.5 * fc_ij * (vr_ij - b_ij * va_ij);
      potential_energy_total = potential_energy_total + potential_energy_partial;
      let fp_j = fp.entry(j_id).or_insert_with(defaultFP);
      fp_j.potential_energy += potential_energy_partial;
    }
  }

  FPInfoBoxed {
    fp,
    potential_energy: potential_energy_total,
    optimization_considered,
    optimization_ignored,
  }
}

impl BoxedWorld {
  pub fn update_verlet_nose_hoover(&mut self, time_step: f64, next_iteration: usize,
                                   desired_temperature: f64, q_effective_mass: f64) {
    let mut next_iteration_atom_container = SimpleAtomContainer::new_fixed_cap(self.num_of_atoms);
    
    for sim_box in self.box_container.read().unwrap().current_boxes().iter() {
      let vel_task = BoxTask::VelocityTask {
        box_container: self.box_container.clone(),
        box_id: sim_box.id(),
        time_step,
        previous_thermostat_epsilon: self.box_container.read().unwrap().current_thermostat_epsilon(),
        current_iteration: self.iteration,
      };
      
      self.tx_task.send(vel_task).unwrap();
    }
    
    unimplemented!("aha");
  }
}