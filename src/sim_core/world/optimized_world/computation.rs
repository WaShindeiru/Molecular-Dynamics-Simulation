use std::collections::HashMap;

use nalgebra::Vector3;

use crate::data::constants::get_constants;
use crate::data::types::get_interaction_type;
use crate::particle::potential::b::g;
use crate::particle::potential::fc::{fc, fc_gradient};
use crate::particle::potential::va::{va, va_gradient};
use crate::particle::potential::vr::{vr, vr_gradient};
use crate::particle::potential::b;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boundary_constraint::periodic::{
  check_position_constraint_periodic, check_position_constraint_periodic_all,
};
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::HalfVelocityResult;
use crate::sim_core::world::cell::{FixedPositionParticle, LinkedCellContainer};
use crate::sim_core::world::computation::FP;
use crate::utils::math::cos_from_vec;

/// Same physics as `crate::sim_core::world::computation::compute_forces_potential`, but
/// operating on plain `{id, position}` pairs (already periodic-adjusted, see
/// `LinkedCellContainer::atoms_for_cell_fixed`/`neighbour_atoms_periodic_fixed_positions`)
/// instead of `Arc`-wrapped particles behind a `dyn ForceComputationOperations` trait object.
/// `particles_i`/`particles_j` are plain slices, so no per-pair `Arc` clone/drop and no
/// per-call `Vec` re-clone of the interaction-partner set, which the trait-object version pays
/// for on every `(i, j)` pair via `particles_j.clone().into_iter()`.
pub fn compute_forces_potential(
  particles_i: &[FixedPositionParticle],
  particles_j: &[FixedPositionParticle],
  container: &LinkedCellContainer,
  fp: &mut Vec<FP>,
  gradients_cache: &mut Vec<Vector3<f64>>,
) -> f64 {
  let zero_fp = FP { force: Vector3::zeros(), potential_energy: 0. };

  for j in particles_j {
    fp[j.id] = zero_fp.clone();
  }

  let mut potential_energy_total: f64 = 0.;

  for i in particles_i {
    let i_id = i.id;
    let i_type = container.particles().get(i_id).unwrap().get_type();

    for j in particles_j {
      let j_id = j.id;
      if j_id == i_id {
        continue;
      }
      let j_type = container.particles().get(j_id).unwrap().get_type();

      let c_ij = get_constants(&get_interaction_type(&i_type, &j_type));

      let r_ij_vec = j.position - i.position;
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

      for k in particles_j {
        let k_id = k.id;
        if k_id == j_id || k_id == i_id {
          continue;
        }
        let k_type = container.particles().get(k_id).unwrap().get_type();

        let c_ik = get_constants(&get_interaction_type(&i_type, &k_type));

        let r_ik_vec = k.position - i.position;
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

      for k in particles_j {
        let k_id = k.id;
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

/// Same as `crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::
/// verlet_noose_hoover_half_velocity_position`, but takes a reference to the (already fully
/// built) particle container plus a slice of ids to process, instead of an owned iterator of
/// `AsRef<Particle>` items. No periodic-adjusted positions are needed here: the half-velocity
/// step only ever looks at a particle relative to itself, never at pairwise distances.
pub fn verlet_noose_hoover_half_velocity_position(
  container: &LinkedCellContainer,
  ids: &[usize],
  time_step: f64,
  previous_thermostat_epsilon: f64,
  current_iteration: usize,
  container_size: &Vector3<f64>,
  edge_condition: EdgeCondition,
) -> HalfVelocityResult {
  let atom_count = ids.len();
  let mut half_velocity_cache: HashMap<usize, Vector3<f64>> = HashMap::with_capacity(atom_count);
  let mut compliance_cache = HashMap::with_capacity(atom_count);
  let mut new_position_cache: HashMap<usize, Vector3<f64>> = HashMap::with_capacity(atom_count);
  let mut thermostat_work_cache: HashMap<usize, f64> = HashMap::with_capacity(atom_count);

  for &i_id in ids {
    let atom_i = container.particles().get(i_id).unwrap();

    let thermostat_difference =
      atom_i.get_acceleration() - previous_thermostat_epsilon * atom_i.get_velocity();
    let half_velocity_i: Vector3<f64> = atom_i.get_velocity() + thermostat_difference * (time_step / 2.0);
    half_velocity_cache.insert(i_id, half_velocity_i);

    let previous_position = atom_i.get_position();
    let next_position: Vector3<f64> = previous_position + half_velocity_i * time_step;

    let (validated_position, compliance) = match edge_condition {
      EdgeCondition::Simple { .. } => unimplemented!("not yet"),
      EdgeCondition::Periodic { .. } => {
        check_position_constraint_periodic(next_position, container_size)
      }
      EdgeCondition::PeriodicAll => check_position_constraint_periodic_all(next_position, container_size),
    };

    let thermostat_work = if current_iteration == 0 || !atom_i.is_atom() {
      0.
    } else {
      let thermostat_force = previous_thermostat_epsilon * atom_i.get_mass() * atom_i.get_velocity();
      let thermostat_path = next_position - previous_position;
      thermostat_force.magnitude() * thermostat_path.magnitude() * cos_from_vec(&thermostat_force, &thermostat_path)
    };

    compliance_cache.insert(i_id, compliance);
    new_position_cache.insert(i_id, validated_position);
    thermostat_work_cache.insert(i_id, thermostat_work);
  }

  HalfVelocityResult {
    half_velocity: half_velocity_cache,
    new_position: new_position_cache,
    compliance: compliance_cache,
    thermostat_work: thermostat_work_cache,
  }
}
