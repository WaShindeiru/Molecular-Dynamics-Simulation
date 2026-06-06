use std::collections::HashSet;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::data::constants::ATOMIC_MASS_FE;
use crate::data::types::AtomType;
use crate::particle::Atom;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::{
  BoxContainerConfig, new_config,
};
use crate::sim_core::world::boxed_world::box_container::get_needed_box_id_periodic;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::force_task_box_container::ForceTaskBoxContainer;

use super::compute_forces_for_boxes_primary_only;

fn fe_atom(id: usize, position: Vector3<f64>) -> Particle {
  Particle::Atom(Atom::new(
    id,
    AtomType::Fe,
    ATOMIC_MASS_FE,
    position,
    Vector3::zeros(),
    Vector3::zeros(),
    Vector3::zeros(),
    0.0,
  ))
}

/// FeFe, world 8³ → 2×2×2 cells, box length 4.
fn test_config() -> BoxContainerConfig {
  let world_size = Vector3::new(8.0, 8.0, 8.0);
  new_config(&[fe_atom(0, Vector3::zeros())], world_size)
}

fn value_container(particles: &[(usize, Vector3<f64>)]) -> BoxContainer<SimulationBox> {
  let mut container = BoxContainer::<SimulationBox>::new_local(test_config());
  for &(id, position) in particles {
    container.add_particle(Arc::new(fe_atom(id, position)));
  }
  container
}

fn option_container(particles: &[(usize, Vector3<f64>)]) -> BoxContainer<Option<Arc<SimulationBox>>> {
  value_container(particles)
    .into_shared()
    .to_force_option_view()
}

fn periodic_edge() -> EdgeCondition {
  EdgeCondition::Periodic {
    trigger_small_subtask_size: 1,
    split: EdgeCondition::DEFAULT_SPLIT,
  }
}

fn build_force_container(
  working: &[(usize, Vector3<f64>)],
  static_halo: &[(usize, Vector3<f64>)],
  primary_box_ids: &[usize],
  force_center_box_ids: &[usize],
) -> ForceTaskBoxContainer {
  let config = test_config();
  let working_boxes = value_container(working);
  let static_view = option_container(static_halo);
  let mut force_view = BoxContainer::merge_force_views(working_boxes, &static_view)
    .expect("working and static halo must not share particle ids");

  let mut needed_boxes: HashSet<usize> =
    get_needed_box_id_periodic(force_center_box_ids, &config)
      .into_iter()
      .collect();
  needed_boxes.extend(primary_box_ids.iter().copied());
  let needed_boxes: Vec<usize> = needed_boxes.into_iter().collect();
  force_view.ensure_boxes_exist(&needed_boxes);

  ForceTaskBoxContainer::new(force_view, periodic_edge())
}

fn working_box_ids_from_particles(particles: &[(usize, Vector3<f64>)]) -> Vec<usize> {
  let config = test_config();
  let mut ids: HashSet<usize> = HashSet::new();
  for &(_, position) in particles {
    ids.insert(config.box_id_for_position(&position));
  }
  let mut ids: Vec<usize> = ids.into_iter().collect();
  ids.sort_unstable();
  ids
}

fn assert_finite_primary_forces(
  forces: &std::collections::HashMap<usize, Vector3<f64>>,
  primary_ids: &HashSet<usize>,
) {
  assert_eq!(forces.len(), primary_ids.len());
  for &id in primary_ids {
    let force = forces.get(&id).expect("primary particle must have a force");
    assert!(force.x.is_finite(), "particle {id} force.x");
    assert!(force.y.is_finite(), "particle {id} force.y");
    assert!(force.z.is_finite(), "particle {id} force.z");
  }
}

fn assert_primary_forces_non_zero_for(
  forces: &std::collections::HashMap<usize, Vector3<f64>>,
  ids: &[usize],
) {
  for &id in ids {
    let force = forces.get(&id).expect("primary particle must have a force");
    assert!(force.norm() > 0.0, "particle {id} should feel non-zero force");
  }
}

/// Primary box (0,0,0): 3 working particles. Neighbours: one with 1 halo atom,
/// one with 2 halo atoms, the rest empty.
#[test]
fn compute_forces_all_three_particles_in_primary_box() {
  let config = test_config();
  let primary_box = config.box_id_for_position(&Vector3::new(0.5, 0.5, 0.5));
  let neighbor_x = config.box_id_for_position(&Vector3::new(4.5, 0.5, 0.5));
  let neighbor_y = config.box_id_for_position(&Vector3::new(0.5, 4.5, 0.5));

  let working = [
    (0, Vector3::new(0.5, 0.5, 0.5)),
    (1, Vector3::new(1.0, 1.0, 1.0)),
    (2, Vector3::new(2.0, 0.5, 1.0)),
  ];
  let static_halo = [
    (100, Vector3::new(4.5, 0.5, 0.5)),
    (101, Vector3::new(0.5, 4.5, 0.5)),
    (102, Vector3::new(1.0, 4.5, 1.0)),
  ];

  let primary_ids: HashSet<usize> = [0, 1, 2].into_iter().collect();
  let force_center_box_ids = working_box_ids_from_particles(&working);
  let force_container =
    build_force_container(&working, &static_halo, &[primary_box], &force_center_box_ids);

  assert_eq!(force_center_box_ids, vec![primary_box]);
  assert!(force_container.view().get_box(neighbor_x).is_some());
  assert!(force_container.view().get_box(neighbor_y).is_some());

  let forces = compute_forces_for_boxes_primary_only(
    &force_container,
    &force_center_box_ids,
    &primary_ids,
  );

  assert_finite_primary_forces(&forces, &primary_ids);
  assert_primary_forces_non_zero_for(&forces, &[0, 1, 2]);
}

/// Primary box holds 2 working particles; the third primary particle sits in +x neighbour.
/// Halo neighbours mix empty / 1 / 2 particles.
#[test]
fn compute_forces_one_primary_particle_in_neighbor_box() {
  let config = test_config();
  let primary_box = config.box_id_for_position(&Vector3::new(0.5, 0.5, 0.5));
  let neighbor_x = config.box_id_for_position(&Vector3::new(4.5, 0.5, 0.5));
  let neighbor_y = config.box_id_for_position(&Vector3::new(0.5, 4.5, 0.5));
  let neighbor_z = config.box_id_for_position(&Vector3::new(0.5, 0.5, 4.5));

  let working = [
    (0, Vector3::new(3.5, 0.5, 0.5)),
    (1, Vector3::new(3.0, 1.0, 1.0)),
    (2, Vector3::new(4.1, 0.5, 0.5)),
  ];
  let static_halo = [
    (100, Vector3::new(0.5, 4.5, 0.5)),
    (101, Vector3::new(0.5, 0.5, 4.5)),
    (102, Vector3::new(1.0, 0.5, 4.5)),
  ];

  let primary_ids: HashSet<usize> = [0, 1, 2].into_iter().collect();
  let force_center_box_ids = working_box_ids_from_particles(&working);
  let force_container =
    build_force_container(&working, &static_halo, &[primary_box], &force_center_box_ids);

  assert_eq!(force_center_box_ids, vec![primary_box, neighbor_x]);
  assert!(force_container.view().get_box(primary_box).is_some());
  assert!(force_container.view().get_box(neighbor_y).is_some());
  assert!(force_container.view().get_box(neighbor_z).is_some());

  let forces = compute_forces_for_boxes_primary_only(
    &force_container,
    &force_center_box_ids,
    &primary_ids,
  );

  assert_finite_primary_forces(&forces, &primary_ids);
  assert_primary_forces_non_zero_for(&forces, &[0, 1, 2]);
}

/// All 3 primary particles left the assigned primary box; it exists empty via ensure_boxes_exist.
/// Particles occupy three different neighbour boxes; remaining neighbours are empty / 1 / 2 halo atoms.
#[test]
fn compute_forces_all_three_particles_in_neighbor_boxes() {
  let config = test_config();
  let primary_box = config.box_id_for_position(&Vector3::new(0.5, 0.5, 0.5));
  let neighbor_x = config.box_id_for_position(&Vector3::new(4.5, 0.5, 0.5));
  let neighbor_xy = config.box_id_for_position(&Vector3::new(4.5, 4.5, 0.5));

  let working = [
    (0, Vector3::new(4.1, 4.1, 0.5)),
    (1, Vector3::new(4.2, 4.15, 0.5)),
    (2, Vector3::new(4.1, 3.95, 0.5)),
  ];
  let static_halo = [
    (100, Vector3::new(0.5, 0.5, 4.5)),
    (101, Vector3::new(4.5, 0.5, 4.5)),
    (102, Vector3::new(4.5, 4.5, 4.5)),
  ];

  let primary_ids: HashSet<usize> = [0, 1, 2].into_iter().collect();
  let force_center_box_ids = working_box_ids_from_particles(&working);
  let force_container =
    build_force_container(&working, &static_halo, &[primary_box], &force_center_box_ids);

  assert_eq!(force_center_box_ids, vec![neighbor_x, neighbor_xy]);
  assert!(force_container.view().get_box(primary_box).is_some());
  assert_eq!(
    force_container
      .view()
      .get_box(primary_box)
      .unwrap()
      .particles()
      .len(),
    0,
    "primary box must exist but stay empty",
  );

  let forces = compute_forces_for_boxes_primary_only(
    &force_container,
    &force_center_box_ids,
    &primary_ids,
  );

  assert_finite_primary_forces(&forces, &primary_ids);
  assert_primary_forces_non_zero_for(&forces, &[0, 1, 2]);
}
