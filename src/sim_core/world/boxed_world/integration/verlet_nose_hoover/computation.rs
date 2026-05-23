use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::periodic::{check_position_constraint_periodic, check_position_constraint_periodic_all, check_position_constraint_periodic_quadratic};
use crate::sim_core::world::boundary_constraint::simple::check_position_constraint_simple;
use crate::sim_core::world::boundary_constraint::{EdgeCondition, ParticleCompliance};
use crate::utils::math::cos_from_vec;
use nalgebra::Vector3;
use std::collections::HashMap;

pub use crate::sim_core::world::computation::{
  FP, FPInfoBoxed, ForceComputationOperations, compute_forces_potential,
};

pub struct HalfVelocityResult {
  pub half_velocity: HashMap<usize, Vector3<f64>>,
  pub new_position: HashMap<usize, Vector3<f64>>,
  pub thermostat_work: HashMap<usize, f64>,
  pub compliance: HashMap<usize, ParticleCompliance>,
}

pub fn verlet_noose_hoover_half_velocity_position<I>(
  previous_atom_container: I,
  time_step: f64,
  previous_thermostat_epsilon: f64,
  atom_count: usize,
  current_iteration: usize,
  container_size: &Vector3<f64>,
  edge_condition: EdgeCondition,
) -> HalfVelocityResult
where
  I: IntoIterator,
  I::Item: AsRef<Particle>,
{
  let mut half_velocity_cache: HashMap<usize, Vector3<f64>> = HashMap::with_capacity(atom_count);
  let mut compliance_cache: HashMap<usize, ParticleCompliance> = HashMap::with_capacity(atom_count);
  let mut new_position_cache: HashMap<usize, Vector3<f64>> = HashMap::with_capacity(atom_count);
  let mut thermostat_work_cache: HashMap<usize, f64> = HashMap::with_capacity(atom_count);

  for temp_i in previous_atom_container.into_iter() {
    let atom_i = temp_i.as_ref();
    let i_id = atom_i.get_id();
    let previous_velocity = atom_i.get_velocity();

    let thermostat_difference =
      atom_i.get_acceleration() - previous_thermostat_epsilon * atom_i.get_velocity();
    let half_velocity_i: Vector3<f64> =
      atom_i.get_velocity() + thermostat_difference * (time_step / 2.0);
    half_velocity_cache.insert(i_id, half_velocity_i);

    let previous_position = atom_i.get_position();
    let next_position: Vector3<f64> = previous_position + half_velocity_i * time_step;

    let (validated_position, compliance) = match edge_condition {
      EdgeCondition::Simple { .. } => {
        unimplemented!("not yet");
        check_position_constraint_simple(next_position, container_size)
      },
      EdgeCondition::Periodic { .. } => check_position_constraint_periodic(
        next_position,
        container_size
      ),
      EdgeCondition::PeriodicAll => check_position_constraint_periodic_all(next_position, container_size),
    };

    let thermostat_work;

    if current_iteration == 0 {
      thermostat_work = 0.;
    } else {
      let thermostat_force =
        previous_thermostat_epsilon * atom_i.get_mass() * atom_i.get_velocity();
      let thermostat_path = next_position - previous_position;
      thermostat_work = thermostat_force.magnitude()
        * thermostat_path.magnitude()
        * cos_from_vec(&thermostat_force, &thermostat_path);
    }

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

