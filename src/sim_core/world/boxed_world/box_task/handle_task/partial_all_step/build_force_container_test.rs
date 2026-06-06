use std::collections::HashSet;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::data::config::correction_param::CorrectionParam;
use crate::data::constants::ATOMIC_MASS_FE;
use crate::data::types::AtomType;
use crate::data::SimulationConfig;
use crate::particle::Atom;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::sim_core::world::boxed_world::box_task::task_manager::TaskManagerConfig;
use crate::sim_core::world::saver::SaveOptions;
use crate::sim_core::world::thermostat::IntegrationAlgorithm;
use crate::sim_core::world::WorldType;

use super::PartialAllStep;

fn fe_atom(id: usize, position: Vector3<f64>) -> Particle {
  let mass = ATOMIC_MASS_FE;
  Particle::Atom(Atom::new(
    id,
    AtomType::Fe,
    mass,
    position,
    Vector3::zeros(),
    Vector3::zeros(),
    Vector3::zeros(),
    0.0,
  ))
}

/// Two adjacent boxes along x (FeFe, world 8³ → 2×2×2 cells, box length 4).
fn two_box_history() -> BoxContainer<Arc<SimulationBox>> {
  let world_size = Vector3::new(8.0, 8.0, 8.0);
  let atoms = vec![
    fe_atom(0, Vector3::new(0.5, 0.5, 0.5)),
    fe_atom(1, Vector3::new(4.5, 0.5, 0.5)),
  ];
  BoxContainer::new(atoms, world_size)
}

fn test_simulation_config(world_size: Vector3<f64>) -> SimulationConfig {
  SimulationConfig::new(
    world_size,
    0.0,
    0.001,
    100,
    1000,
    SaveOptions::default(),
    IntegrationAlgorithm::VelocityVerlet,
    WorldType::BoxedWorld {
      task_manager_config: TaskManagerConfig {
        debug: false,
        task_worker_multiplier: 4.0,
      },
    },
    EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    },
    CorrectionParam::default(),
  )
}

#[test]
fn merge_contains_all_particles_with_correct_box_id_mapping() {
  let history = Arc::new(two_box_history());
  let config = test_simulation_config(history.config().world_size);

  println!("{}", history.config().box_count_dim);

  let box_a = history
    .config()
    .box_id_for_position(&Vector3::new(0.5, 0.5, 0.5));

  let step = PartialAllStep::new(&[box_a], history.clone(), config, 0.0, 0);
  let force_container = step.build_force_container_for_test();
  let view = force_container.view();

  let expected_ids: HashSet<usize> = history
    .all_particles()
    .map(|particle| particle.get_id())
    .collect();
  let merged_ids: HashSet<usize> = view.box_id_cache().keys().copied().collect();
  assert_eq!(merged_ids, expected_ids);

  for (&particle_id, &mapped_box_id) in view.box_id_cache() {
    let sim_box = view
      .get_box(mapped_box_id)
      .expect("mapped box must exist in force view");
    let particle = sim_box
      .particles()
      .get(&particle_id)
      .expect("particle must live in its mapped box");
    let expected_box_id = view
      .config()
      .box_id_for_position(particle.get_position());
    assert_eq!(mapped_box_id, expected_box_id);
  }
}

#[test]
fn working_particles_come_from_working_not_history_v1() {
  let history = Arc::new(two_box_history());
  let config = test_simulation_config(history.config().world_size);
  let box_a = history
    .config()
    .box_id_for_position(&Vector3::new(0.5, 0.5, 0.5));

  let work_particle_id = 0;
  let halo_particle_id = 1;

  let history_position = Vector3::new(0.5, 0.5, 0.5);
  let working_position = Vector3::new(1.2, 0.5, 0.5);
  let static_halo_position = Vector3::new(4.5, 0.5, 0.5);

  let history_particle = history.get_particle(work_particle_id);
  assert_eq!(history_particle.get_position(), &history_position);

  let mut step = PartialAllStep::new(&[box_a], history, config, 0.0, 0);
  step.set_working_position_for_test(0, working_position);

  let force_container = step.build_force_container_for_test();
  let view = force_container.view();

  let primary_box = view
    .get_box(box_a)
    .expect("primary box must exist");
  let work_particle_box_id = view
    .particle_box_id(work_particle_id);
  assert_eq!(work_particle_box_id, box_a);
  let primary = primary_box
    .particles()
    .get(&0)
    .expect("working primary must be in force view");
  assert_eq!(primary.get_position(), &working_position);
  assert_ne!(primary.get_position(), &history_position);

  let static_halo_box = view
    .config()
    .box_id_for_position(&static_halo_position);
  let halo_particle_box_id = view
    .particle_box_id(halo_particle_id);
  assert_eq!(static_halo_box, halo_particle_box_id);
  let halo_box = view
    .get_box(static_halo_box)
    .expect("static halo box must exist");
  let halo = halo_box
    .particles()
    .get(&1)
    .expect("static halo particle must be merged from history");
  assert_eq!(halo.get_position(), &static_halo_position);
}

/// Two adjacent boxes along x (CC box length 2, world 4³).
fn two_box_history_v2() -> BoxContainer<Arc<SimulationBox>> {
  let world_size = Vector3::new(8.0, 8.0, 8.0);
  let atoms = vec![
    fe_atom(0, Vector3::new(4.6, 0.5, 0.5)),
    fe_atom(1, Vector3::new(1.8, 0.5, 0.5)),
  ];
  BoxContainer::new(atoms, world_size)
}

#[test]
fn working_particles_come_from_working_not_history_v2() {
  let history_position = Vector3::new(4.6, 0.5, 0.5);
  let working_position = Vector3::new(1.2, 0.5, 0.5);
  let static_halo_position = Vector3::new(1.8, 0.5, 0.5);

  let history = Arc::new(two_box_history_v2());
  let config = test_simulation_config(history.config().world_size);
  let box_a = history
    .config()
    .box_id_for_position(&history_position);

  let box_b = history
    .config()
    .box_id_for_position(&static_halo_position);

  let work_particle_id = 0;
  let halo_particle_id = 1;

  let history_particle = history.get_particle(work_particle_id);
  assert_eq!(history_particle.get_position(), &history_position);

  let mut step = PartialAllStep::new(&[box_a], history, config, 0.0, 0);
  step.set_working_position_for_test(work_particle_id, working_position);

  let force_container = step.build_force_container_for_test();
  let view = force_container.view();

  let primary_box = view
    .get_box(box_a)
    .expect("primary box must exist");

  let work_particle_box_id = view
    .particle_box_id(work_particle_id);
  assert_eq!(work_particle_box_id, box_b);

  let work_box_id = view
    .particle_box_id(work_particle_id);

  let work_box = view
    .get_box(work_box_id)
    .unwrap();
  
  let mut visited = false;
  for (id, particle) in work_box.particles() {
    if *id == work_particle_id {
      visited = true;
      assert_eq!(particle.get_position(), &working_position);
    }
  }
  assert!(visited);

  let static_halo_box = view
    .config()
    .box_id_for_position(&static_halo_position);
  let halo_particle_box_id = view
    .particle_box_id(halo_particle_id);
  assert_eq!(static_halo_box, halo_particle_box_id);

  let halo_box = view
    .get_box(static_halo_box)
    .expect("static halo box must exist");
  let halo = halo_box
    .particles()
    .get(&1)
    .expect("static halo particle must be merged from history");
  assert_eq!(halo.get_position(), &static_halo_position);
}
