use std::collections::HashSet;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::data::constants::ATOMIC_MASS_FE;
use crate::data::types::AtomType;
use crate::particle::Atom;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::{
  BoxContainerConfig, new_config,
};
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;

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

fn count_particles_in_view(view: &BoxContainer<Option<Arc<SimulationBox>>>) -> usize {
  view
    .simulation_boxes
    .iter()
    .filter_map(|cell| cell.as_ref())
    .map(|sim_box| sim_box.particles().len())
    .sum()
}

fn particle_position_in_value(
  container: &BoxContainer<SimulationBox>,
  particle_id: usize,
) -> Vector3<f64> {
  let box_id = container.particle_box_id(particle_id);
  *container
    .get_box(box_id)
    .particles()
    .get(&particle_id)
    .expect("particle must exist in value container")
    .get_position()
}

fn particle_position_in_option(
  container: &BoxContainer<Option<Arc<SimulationBox>>>,
  particle_id: usize,
) -> Vector3<f64> {
  let box_id = container.particle_box_id(particle_id);
  *container
    .get_box(box_id)
    .expect("box must exist in option container")
    .particles()
    .get(&particle_id)
    .expect("particle must exist in option container")
    .get_position()
}

fn particle_position_in_merged(
  merged: &BoxContainer<Option<Arc<SimulationBox>>>,
  particle_id: usize,
) -> Vector3<f64> {
  let box_id = merged.particle_box_id(particle_id);
  *merged
    .get_box(box_id)
    .expect("mapped box must exist in merged view")
    .particles()
    .get(&particle_id)
    .expect("particle must live in its mapped box")
    .get_position()
}

fn assert_particles_match(
  merged: &BoxContainer<Option<Arc<SimulationBox>>>,
  expected: &[(usize, Vector3<f64>)],
) {
  assert_eq!(merged.box_id_cache().len(), expected.len());
  assert_eq!(count_particles_in_view(merged), expected.len());

  let expected_ids: HashSet<usize> = expected.iter().map(|(id, _)| *id).collect();
  let merged_ids: HashSet<usize> = merged.box_id_cache().keys().copied().collect();
  assert_eq!(merged_ids, expected_ids);

  for &(particle_id, expected_position) in expected {
    let mapped_box_id = merged.particle_box_id(particle_id);
    let expected_box_id = merged
      .config()
      .box_id_for_position(&expected_position);
    assert_eq!(mapped_box_id, expected_box_id);
    assert_eq!(
      particle_position_in_merged(merged, particle_id),
      expected_position,
    );
  }
}

#[test]
fn merge_combines_disjoint_working_and_static_particles() {
  let primary_position = Vector3::new(0.5, 0.5, 0.5);
  let halo_position = Vector3::new(4.5, 0.5, 0.5);

  let working = value_container(&[(0, primary_position)]);
  let static_halo = option_container(&[(1, halo_position)]);
  let working_source_position = particle_position_in_value(&working, 0);
  let static_source_position = particle_position_in_option(&static_halo, 1);

  let merged = BoxContainer::merge_force_views(working, &static_halo)
    .expect("disjoint sources should merge");

  assert_particles_match(&merged, &[(0, primary_position), (1, halo_position)]);
  assert_eq!(
    particle_position_in_merged(&merged, 0),
    working_source_position,
  );
  assert_eq!(
    particle_position_in_merged(&merged, 1),
    static_source_position,
  );
}

#[test]
fn merge_preserves_working_positions_when_static_fills_neighbor_box() {
  let primary_position = Vector3::new(0.5, 0.5, 0.5);
  let halo_position = Vector3::new(1.2, 0.5, 0.5);
  let box_a = test_config().box_id_for_position(&primary_position);

  let working = value_container(&[(0, primary_position)]);
  let static_halo = option_container(&[(1, halo_position)]);
  let working_source_position = particle_position_in_value(&working, 0);
  let static_source_position = particle_position_in_option(&static_halo, 1);

  let merged = BoxContainer::merge_force_views(working, &static_halo)
    .expect("disjoint sources should merge");

  assert_particles_match(&merged, &[(0, primary_position), (1, halo_position)]);
  assert_eq!(
    particle_position_in_merged(&merged, 0),
    working_source_position,
  );
  assert_eq!(
    particle_position_in_merged(&merged, 1),
    static_source_position,
  );

  let primary_box = merged
    .get_box(box_a)
    .expect("primary box must exist");
  assert_eq!(
    *primary_box.particles().get(&0).unwrap().get_position(),
    primary_position,
  );
  assert_eq!(
    *primary_box.particles().get(&1).unwrap().get_position(),
    halo_position,
  );
}

#[test]
fn merge_returns_none_when_particle_ids_overlap() {
  let primary_position = Vector3::new(0.5, 0.5, 0.5);

  let working = value_container(&[(0, primary_position)]);
  let static_halo = option_container(&[(0, primary_position)]);

  assert!(BoxContainer::merge_force_views(working, &static_halo).is_none());
}
