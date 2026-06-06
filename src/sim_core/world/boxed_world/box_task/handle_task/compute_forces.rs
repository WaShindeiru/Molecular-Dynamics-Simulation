use std::collections::{HashMap, HashSet};

use crate::sim_core::world::boxed_world::box_task::force_task_box_container::ForceTaskBoxContainer;
use crate::sim_core::world::computation::{ForceComputationOperations, compute_forces_potential};

/// Computes pair forces for each primary box in `primary_box_ids` using `force_container`.
/// Returns accumulated force per particle id (all ids appearing in `particles_j` for those boxes).
pub fn compute_forces_for_boxes(
  force_container: &ForceTaskBoxContainer,
  primary_box_ids: &[usize],
) -> (
  HashMap<usize, crate::sim_core::world::boxed_world::box_task::ForceTaskParticleData>,
  f64,
) {
  let mut particles: HashMap<usize, crate::sim_core::world::boxed_world::box_task::ForceTaskParticleData> =
    HashMap::new();
  let mut potential_energy_total = 0.0f64;

  for &box_id in primary_box_ids {
    let particles_i: Vec<Box<dyn ForceComputationOperations>> = force_container
      .atoms_for_box(box_id)
      .expect("Primary box must be present in ForceTaskBoxContainer")
      .collect();

    let mut particles_j: Vec<Box<dyn ForceComputationOperations>> =
      force_container.neighbour_atoms_periodic(box_id).collect();
    particles_j.extend(force_container.atoms_for_box(box_id).unwrap());

    let info = compute_forces_potential(&particles_i, &particles_j);
    potential_energy_total += info.potential_energy;

    for particle in particles_j.iter() {
      let id = particle.get_id();
      let fp = info.fp.get(&id).unwrap();
      let particle_box_id = force_container.view().particle_box_id(id);
      particles
        .entry(id)
        .and_modify(|entry| {
          entry.force += fp.force;
          entry.potential_energy += fp.potential_energy;
        })
        .or_insert(crate::sim_core::world::boxed_world::box_task::ForceTaskParticleData {
          box_id: particle_box_id,
          force: fp.force,
          potential_energy: fp.potential_energy,
        });
    }
  }

  (particles, potential_energy_total)
}

/// Like `compute_forces_for_boxes` but only returns forces for ids in `primary_ids`.
pub fn compute_forces_for_boxes_primary_only(
  force_container: &ForceTaskBoxContainer,
  primary_box_ids: &[usize],
  primary_ids: &HashSet<usize>,
) -> HashMap<usize, nalgebra::Vector3<f64>> {
  let (all_forces, _) = compute_forces_for_boxes(force_container, primary_box_ids);
  all_forces
    .into_iter()
    .filter(|(id, _)| primary_ids.contains(id))
    .map(|(id, data)| (id, data.force))
    .collect()
}

#[cfg(test)]
mod compute_forces_primary_only_test;
