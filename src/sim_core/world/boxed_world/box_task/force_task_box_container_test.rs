use crate::data::InteractionType::FeC;
  use crate::data::types::AtomType::{C, Fe};
  use crate::data::InteractionType::FeFe;
  use crate::particle::SafeAtomFactory;
  use crate::sim_core::world::computation::compute_forces_potential;
  use super::*;

  fn test_box_container_config() -> BoxContainerConfig {
    BoxContainerConfig {
      box_type: FeC,
      box_length: Vector3::new(10., 10., 10.),
      box_count: 100,
      box_count_dim: Vector3::new(10, 10, 10),
      world_size: Vector3::new(100., 100., 100.),
    }
  }

  #[test]
  fn test_force_task_box_container_from_config_and_needed_box_ids() {
    let container_config = test_box_container_config();

    let mut needed_box_ids: Vec<usize> = Vec::new();

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..8 {
          needed_box_ids.push(get_id_simulation_box(
            &Vector3::new(x, y, z),
            &container_config.box_count_dim,
          ));
        }
      }
    }

    let box_container = BoxContainer::new_local(container_config).into_shared();
    let view = box_container.view_select_boxes(&needed_box_ids);
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    });

    let mut some_box_count = 0;
    let mut none_box_count = 0;

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let box_id =
            get_id_simulation_box(&Vector3::new(x, y, z), &container_config.box_count_dim);

          if force_container.view().get_box(box_id).is_some() {
            some_box_count += 1;
          } else {
            none_box_count += 1;
          }
        }
      }
    }

    assert_eq!(some_box_count, needed_box_ids.len());
    assert_eq!(
      none_box_count,
      container_config.box_count_dim.x * container_config.box_count_dim.y
        * container_config.box_count_dim.z
        - needed_box_ids.len()
    );
  }

  #[test]
  fn test_force_task_box_container_with_one_particle_in_each_box_center_atoms_for_box_normal_placement() {
    let container_config = test_box_container_config();

    let mut needed_box_ids: Vec<usize> = Vec::new();

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..8 {
          needed_box_ids.push(get_id_simulation_box(
            &Vector3::new(x, y, z),
            &container_config.box_count_dim,
          ));
        }
      }
    }

    let mut box_container = BoxContainer::new_local(container_config);
    let atom_factory = SafeAtomFactory::new(1.0, container_config.world_size.z);

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let position = Vector3::new(
            (x as f64 + 0.5) * container_config.box_length.x,
            (y as f64 + 0.5) * container_config.box_length.y,
            (z as f64 + 0.5) * container_config.box_length.z,
          );
          let particle = atom_factory.get_atom(
            C,
            position,
            Vector3::zeros(),
          );

          box_container.add_particle(Arc::new(particle));
        }
      }
    }

    let view = box_container
      .into_shared()
      .view_select_boxes(&needed_box_ids);
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    });

    let mut some_box_count = 0;
    let mut none_box_count = 0;
    let mut particle_count = 0;

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let box_id =
            get_id_simulation_box(&Vector3::new(x, y, z), &container_config.box_count_dim);

          if let Some(sim_box) = force_container.view().get_box(box_id) {
            some_box_count += 1;
            particle_count += sim_box.len();
          } else {
            none_box_count += 1;
          }
        }
      }
    }

    assert_eq!(some_box_count, needed_box_ids.len());
    assert_eq!(
      none_box_count,
      container_config.box_count_dim.x * container_config.box_count_dim.y
        * container_config.box_count_dim.z
        - needed_box_ids.len()
    );
    assert_eq!(particle_count, needed_box_ids.len());

    let box_id = get_id_simulation_box(
      &Vector3::new(3, 3, 3),
      &container_config.box_count_dim,
    );
    let particles: Vec<Box<dyn ForceComputationOperations>> = force_container
      .atoms_for_box(box_id)
      .expect("Box (3, 3, 3) should be present in ForceTaskBoxContainer")
      .collect();

    assert_eq!(particles.len(), 1);

    let particle = particles.first().unwrap();
    let expected_position = Vector3::new(35.0, 35.0, 35.0);
    assert_eq!(particle.get_position(), expected_position);

    let proxy = particle
      .as_any()
      .downcast_ref::<ParticlePositionProxy>()
      .expect("atoms_for_box should return particle position proxies");
    assert_eq!(
      *proxy.particle_placement(),
      ParticlePlacement {
        x: AxisPlacement::Normal,
        y: AxisPlacement::Normal,
        z: AxisPlacement::Normal,
      }
    );
  }

  #[test]
  fn test_force_task_box_container_with_one_particle_in_each_box_center_neighbour_atoms_periodic_normal_placement() {
    let container_config = test_box_container_config();

    let mut needed_box_ids: Vec<usize> = Vec::new();

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..8 {
          needed_box_ids.push(get_id_simulation_box(
            &Vector3::new(x, y, z),
            &container_config.box_count_dim,
          ));
        }
      }
    }

    let mut box_container = BoxContainer::new_local(container_config);
    let atom_factory = SafeAtomFactory::new(1.0, container_config.world_size.z);

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let position = Vector3::new(
            (x as f64 + 0.5) * container_config.box_length.x,
            (y as f64 + 0.5) * container_config.box_length.y,
            (z as f64 + 0.5) * container_config.box_length.z,
          );
          let particle = atom_factory.get_atom(C, position, Vector3::zeros());

          box_container.add_particle(Arc::new(particle));
        }
      }
    }

    let view = box_container
      .into_shared()
      .view_select_boxes(&needed_box_ids);
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    });
    let box_id = get_id_simulation_box(
      &Vector3::new(3, 3, 3),
      &container_config.box_count_dim,
    );

    let particles: Vec<Box<dyn ForceComputationOperations>> =
      force_container.neighbour_atoms_periodic(box_id).collect();

    assert_eq!(particles.len(), 26);

    let mut expected_positions: Vec<Vector3<f64>> = Vec::new();
    for x in 2..=4 {
      for y in 2..=4 {
        for z in 2..=4 {
          if x == 3 && y == 3 && z == 3 {
            continue;
          }

          expected_positions.push(Vector3::new(
            (x as f64 + 0.5) * container_config.box_length.x,
            (y as f64 + 0.5) * container_config.box_length.y,
            (z as f64 + 0.5) * container_config.box_length.z,
          ));
        }
      }
    }

    for particle in particles {
      let proxy = particle
        .as_any()
        .downcast_ref::<ParticlePositionProxy>()
        .expect("neighbour_atoms_periodic should return particle position proxies");
      assert_eq!(
        *proxy.particle_placement(),
        ParticlePlacement {
          x: AxisPlacement::Normal,
          y: AxisPlacement::Normal,
          z: AxisPlacement::Normal,
        }
      );

      let position = particle.get_position();
      let expected_index = expected_positions
        .iter()
        .position(|expected_position| *expected_position == position)
        .expect("Neighbour particle should have expected center position");
      expected_positions.remove(expected_index);
    }

    assert!(expected_positions.is_empty());
  }

  #[test]
  fn test_force_task_box_container_with_one_particle_in_each_box_center_neighbour_atoms_periodic_left_edge_placement() {
    let container_config = test_box_container_config();

    let mut needed_box_ids: Vec<usize> = Vec::new();

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..8 {
          needed_box_ids.push(get_id_simulation_box(
            &Vector3::new(x, y, z),
            &container_config.box_count_dim,
          ));
        }
      }
    }

    let mut box_container = BoxContainer::new_local(container_config);
    let atom_factory = SafeAtomFactory::new(1.0, container_config.world_size.z);

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let position = Vector3::new(
            (x as f64 + 0.5) * container_config.box_length.x,
            (y as f64 + 0.5) * container_config.box_length.y,
            (z as f64 + 0.5) * container_config.box_length.z,
          );
          let particle = atom_factory.get_atom(C, position, Vector3::zeros());

          box_container.add_particle(Arc::new(particle));
        }
      }
    }

    let view = box_container
      .into_shared()
      .view_select_boxes(&needed_box_ids);
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    });
    let box_id = get_id_simulation_box(
      &Vector3::new(0, 3, 3),
      &container_config.box_count_dim,
    );

    let particles: Vec<Box<dyn ForceComputationOperations>> =
      force_container.neighbour_atoms_periodic(box_id).collect();

    assert_eq!(particles.len(), 26);

    let normal_placement = ParticlePlacement {
      x: AxisPlacement::Normal,
      y: AxisPlacement::Normal,
      z: AxisPlacement::Normal,
    };
    let left_x_placement = ParticlePlacement {
      x: AxisPlacement::Left,
      y: AxisPlacement::Normal,
      z: AxisPlacement::Normal,
    };
    let mut expected_particles: Vec<(Vector3<f64>, ParticlePlacement)> = vec![
      (Vector3::new(95.0, 25.0, 25.0), normal_placement),
      (Vector3::new(95.0, 25.0, 35.0), normal_placement),
      (Vector3::new(95.0, 25.0, 45.0), normal_placement),
      (Vector3::new(95.0, 35.0, 25.0), normal_placement),
      (Vector3::new(95.0, 35.0, 35.0), normal_placement),
      (Vector3::new(95.0, 35.0, 45.0), normal_placement),
      (Vector3::new(95.0, 45.0, 25.0), normal_placement),
      (Vector3::new(95.0, 45.0, 35.0), normal_placement),
      (Vector3::new(95.0, 45.0, 45.0), normal_placement),
      (Vector3::new(105.0, 25.0, 25.0), left_x_placement),
      (Vector3::new(105.0, 25.0, 35.0), left_x_placement),
      (Vector3::new(105.0, 25.0, 45.0), left_x_placement),
      (Vector3::new(105.0, 35.0, 25.0), left_x_placement),
      (Vector3::new(105.0, 35.0, 45.0), left_x_placement),
      (Vector3::new(105.0, 45.0, 25.0), left_x_placement),
      (Vector3::new(105.0, 45.0, 35.0), left_x_placement),
      (Vector3::new(105.0, 45.0, 45.0), left_x_placement),
      (Vector3::new(115.0, 25.0, 25.0), left_x_placement),
      (Vector3::new(115.0, 25.0, 35.0), left_x_placement),
      (Vector3::new(115.0, 25.0, 45.0), left_x_placement),
      (Vector3::new(115.0, 35.0, 25.0), left_x_placement),
      (Vector3::new(115.0, 35.0, 35.0), left_x_placement),
      (Vector3::new(115.0, 35.0, 45.0), left_x_placement),
      (Vector3::new(115.0, 45.0, 25.0), left_x_placement),
      (Vector3::new(115.0, 45.0, 35.0), left_x_placement),
      (Vector3::new(115.0, 45.0, 45.0), left_x_placement),
    ];

    for particle in particles {
      let proxy = particle
        .as_any()
        .downcast_ref::<ParticlePositionProxy>()
        .expect("neighbour_atoms_periodic should return particle position proxies");
      let actual_position = particle.get_position();
      let actual_placement = *proxy.particle_placement();
      let expected_index = expected_particles
        .iter()
        .position(|(expected_position, expected_placement)| {
          *expected_position == actual_position && *expected_placement == actual_placement
        })
        .expect("Neighbour particle should have expected center position and placement");

      expected_particles.remove(expected_index);
    }

    assert!(expected_particles.is_empty());
  }

  #[test]
  fn test_force_task_box_container_neighbour_atoms_periodic_corner_left_edge_placement() {
    let container_config = test_box_container_config();

    // Home box at grid (0, y-1, z-1) with box_count_dim = (10, 10, 10) -> (0, 9, 9).
    let home_box_id = get_id_simulation_box(
      &Vector3::new(0, 9, 9),
      &container_config.box_count_dim,
    );
    let needed_box_ids = get_needed_box_id_periodic(&vec![home_box_id], &container_config);

    let mut box_container = BoxContainer::new_local(container_config);
    let atom_factory = SafeAtomFactory::new(1.0, container_config.world_size.z);

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let position = Vector3::new(
            (x as f64 + 0.5) * container_config.box_length.x,
            (y as f64 + 0.5) * container_config.box_length.y,
            (z as f64 + 0.5) * container_config.box_length.z,
          );
          let particle = atom_factory.get_atom(C, position, Vector3::zeros());
          box_container.add_particle(Arc::new(particle));
        }
      }
    }

    let view = box_container
      .into_shared()
      .view_select_boxes(&needed_box_ids);
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    });

    let particles: Vec<Box<dyn ForceComputationOperations>> =
      force_container.neighbour_atoms_periodic(home_box_id).collect();

    assert_eq!(particles.len(), 26);
    assert_eq!(needed_box_ids.len(), 27);

    let normal_placement = ParticlePlacement {
      x: AxisPlacement::Normal,
      y: AxisPlacement::Normal,
      z: AxisPlacement::Normal,
    };
    let left_x_placement = ParticlePlacement {
      x: AxisPlacement::Left,
      y: AxisPlacement::Normal,
      z: AxisPlacement::Normal,
    };
    let left_y_placement = ParticlePlacement {
      x: AxisPlacement::Normal,
      y: AxisPlacement::Left,
      z: AxisPlacement::Normal,
    };
    let left_xy_placement = ParticlePlacement {
      x: AxisPlacement::Left,
      y: AxisPlacement::Left,
      z: AxisPlacement::Normal,
    };

    let mut expected_particles: Vec<(Vector3<f64>, ParticlePlacement)> = vec![
      // x-offset -1 (wrapped to x=9)
      (Vector3::new(95.0, 85.0, 85.0), normal_placement),
      (Vector3::new(95.0, 85.0, 95.0), normal_placement),
      (Vector3::new(95.0, 85.0, 5.0), normal_placement),
      (Vector3::new(95.0, 95.0, 85.0), normal_placement),
      (Vector3::new(95.0, 95.0, 95.0), normal_placement),
      (Vector3::new(95.0, 95.0, 5.0), normal_placement),
      (Vector3::new(95.0, 105.0, 85.0), left_y_placement),
      (Vector3::new(95.0, 105.0, 95.0), left_y_placement),
      (Vector3::new(95.0, 105.0, 5.0), left_y_placement),
      // x-offset 0 (same x column as home)
      (Vector3::new(105.0, 85.0, 85.0), left_x_placement),
      (Vector3::new(105.0, 85.0, 95.0), left_x_placement),
      (Vector3::new(105.0, 85.0, 5.0), left_x_placement),
      (Vector3::new(105.0, 95.0, 85.0), left_x_placement),
      (Vector3::new(105.0, 95.0, 5.0), left_x_placement),
      (Vector3::new(105.0, 105.0, 85.0), left_xy_placement),
      (Vector3::new(105.0, 105.0, 95.0), left_xy_placement),
      (Vector3::new(105.0, 105.0, 5.0), left_xy_placement),
      // x-offset +1
      (Vector3::new(115.0, 85.0, 85.0), left_x_placement),
      (Vector3::new(115.0, 85.0, 95.0), left_x_placement),
      (Vector3::new(115.0, 85.0, 5.0), left_x_placement),
      (Vector3::new(115.0, 95.0, 85.0), left_x_placement),
      (Vector3::new(115.0, 95.0, 95.0), left_x_placement),
      (Vector3::new(115.0, 95.0, 5.0), left_x_placement),
      (Vector3::new(115.0, 105.0, 85.0), left_xy_placement),
      (Vector3::new(115.0, 105.0, 95.0), left_xy_placement),
      (Vector3::new(115.0, 105.0, 5.0), left_xy_placement),
    ];

    for particle in particles {
      let proxy = particle
        .as_any()
        .downcast_ref::<ParticlePositionProxy>()
        .expect("neighbour_atoms_periodic should return particle position proxies");
      let actual_position = particle.get_position();
      let actual_placement = *proxy.particle_placement();
      let expected_index = expected_particles
        .iter()
        .position(|(expected_position, expected_placement)| {
          *expected_position == actual_position && *expected_placement == actual_placement
        })
        .unwrap_or_else(|| {
          panic!(
            "Neighbour particle not matched: position ({}, {}, {}), placement {:?}",
            actual_position.x,
            actual_position.y,
            actual_position.z,
            actual_placement,
          );
        });

      expected_particles.remove(expected_index);
    }

    assert!(expected_particles.is_empty());
  }

  #[test]
  fn test_force_task_box_container_with_one_particle_in_each_box_center_atoms_for_box_left_edge_placement() {
    let container_config = test_box_container_config();

    let mut needed_box_ids: Vec<usize> = Vec::new();

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..8 {
          needed_box_ids.push(get_id_simulation_box(
            &Vector3::new(x, y, z),
            &container_config.box_count_dim,
          ));
        }
      }
    }

    let mut box_container = BoxContainer::new_local(container_config);
    let atom_factory = SafeAtomFactory::new(1.0, container_config.world_size.z);

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let position = Vector3::new(
            (x as f64 + 0.5) * container_config.box_length.x,
            (y as f64 + 0.5) * container_config.box_length.y,
            (z as f64 + 0.5) * container_config.box_length.z,
          );
          let particle = atom_factory.get_atom(
            C,
            position,
            Vector3::zeros(),
          );

          box_container.add_particle(Arc::new(particle));
        }
      }
    }

    let view = box_container
      .into_shared()
      .view_select_boxes(&needed_box_ids);
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    });

    let mut some_box_count = 0;
    let mut none_box_count = 0;
    let mut particle_count = 0;

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let box_id =
            get_id_simulation_box(&Vector3::new(x, y, z), &container_config.box_count_dim);

          if let Some(sim_box) = force_container.view().get_box(box_id) {
            some_box_count += 1;
            particle_count += sim_box.len();
          } else {
            none_box_count += 1;
          }
        }
      }
    }

    assert_eq!(some_box_count, needed_box_ids.len());
    assert_eq!(
      none_box_count,
      container_config.box_count_dim.x * container_config.box_count_dim.y
        * container_config.box_count_dim.z
        - needed_box_ids.len()
    );
    assert_eq!(particle_count, needed_box_ids.len());

    let box_id = get_id_simulation_box(
      &Vector3::new(0, 5, 5),
      &container_config.box_count_dim,
    );
    let particles: Vec<Box<dyn ForceComputationOperations>> = force_container
      .atoms_for_box(box_id)
      .expect("Box (3, 3, 3) should be present in ForceTaskBoxContainer")
      .collect();

    assert_eq!(particles.len(), 1);

    let particle = particles.first().unwrap();
    let expected_position = Vector3::new(105., 55.0, 55.0);
    assert_eq!(particle.get_position(), expected_position);

    let proxy = particle
      .as_any()
      .downcast_ref::<ParticlePositionProxy>()
      .expect("atoms_for_box should return particle position proxies");
    assert_eq!(
      *proxy.particle_placement(),
      ParticlePlacement {
        x: AxisPlacement::Left,
        y: AxisPlacement::Normal,
        z: AxisPlacement::Normal,
      }
    );

    let box_id = get_id_simulation_box(
      &Vector3::new(5, 0, 5),
      &container_config.box_count_dim,
    );
    let particles: Vec<Box<dyn ForceComputationOperations>> = force_container
      .atoms_for_box(box_id)
      .expect("Box (3, 3, 3) should be present in ForceTaskBoxContainer")
      .collect();

    assert_eq!(particles.len(), 1);

    let particle = particles.first().unwrap();
    let expected_position = Vector3::new(55., 105.0, 55.0);
    assert_eq!(particle.get_position(), expected_position);

    let proxy = particle
      .as_any()
      .downcast_ref::<ParticlePositionProxy>()
      .expect("atoms_for_box should return particle position proxies");
    assert_eq!(
      *proxy.particle_placement(),
      ParticlePlacement {
        x: AxisPlacement::Normal,
        y: AxisPlacement::Left,
        z: AxisPlacement::Normal,
      }
    );

    let box_id = get_id_simulation_box(
      &Vector3::new(0, 0, 5),
      &container_config.box_count_dim,
    );
    let particles: Vec<Box<dyn ForceComputationOperations>> = force_container
      .atoms_for_box(box_id)
      .expect("Box (3, 3, 3) should be present in ForceTaskBoxContainer")
      .collect();

    assert_eq!(particles.len(), 1);

    let particle = particles.first().unwrap();
    let expected_position = Vector3::new(105., 105.0, 55.0);
    assert_eq!(particle.get_position(), expected_position);

    let proxy = particle
      .as_any()
      .downcast_ref::<ParticlePositionProxy>()
      .expect("atoms_for_box should return particle position proxies");
    assert_eq!(
      *proxy.particle_placement(),
      ParticlePlacement {
        x: AxisPlacement::Left,
        y: AxisPlacement::Left,
        z: AxisPlacement::Normal,
      }
    );
  }

  #[test]
  fn test_force_task_box_container_with_one_particle_in_each_box_center_atoms_for_box_right_edge_placement() {
    let container_config = test_box_container_config();

    let mut needed_box_ids: Vec<usize> = Vec::new();

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..8 {
          needed_box_ids.push(get_id_simulation_box(
            &Vector3::new(x, y, z),
            &container_config.box_count_dim,
          ));
        }
      }
    }

    let mut box_container = BoxContainer::new_local(container_config);
    let atom_factory = SafeAtomFactory::new(1.0, container_config.world_size.z);

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let position = Vector3::new(
            (x as f64 + 0.5) * container_config.box_length.x,
            (y as f64 + 0.5) * container_config.box_length.y,
            (z as f64 + 0.5) * container_config.box_length.z,
          );
          let particle = atom_factory.get_atom(
            C,
            position,
            Vector3::zeros(),
          );

          box_container.add_particle(Arc::new(particle));
        }
      }
    }

    let view = box_container
      .into_shared()
      .view_select_boxes(&needed_box_ids);
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic {
      trigger_small_subtask_size: 1,
      split: EdgeCondition::DEFAULT_SPLIT,
    });

    let mut some_box_count = 0;
    let mut none_box_count = 0;
    let mut particle_count = 0;

    for x in 0..container_config.box_count_dim.x {
      for y in 0..container_config.box_count_dim.y {
        for z in 0..container_config.box_count_dim.z {
          let box_id =
            get_id_simulation_box(&Vector3::new(x, y, z), &container_config.box_count_dim);

          if let Some(sim_box) = force_container.view().get_box(box_id) {
            some_box_count += 1;
            particle_count += sim_box.len();
          } else {
            none_box_count += 1;
          }
        }
      }
    }

    assert_eq!(some_box_count, needed_box_ids.len());
    assert_eq!(
      none_box_count,
      container_config.box_count_dim.x * container_config.box_count_dim.y
        * container_config.box_count_dim.z
        - needed_box_ids.len()
    );
    assert_eq!(particle_count, needed_box_ids.len());

    let box_id = get_id_simulation_box(
      &Vector3::new(0, 5, 5),
      &container_config.box_count_dim,
    );
    let particles: Vec<Box<dyn ForceComputationOperations>> = force_container
      .atoms_for_box(box_id)
      .expect("Box (3, 3, 3) should be present in ForceTaskBoxContainer")
      .collect();

    assert_eq!(particles.len(), 1);

    let box_id = get_id_simulation_box(
      &Vector3::new(5, 9, 5),
      &container_config.box_count_dim,
    );
    let particles: Vec<Box<dyn ForceComputationOperations>> = force_container
      .atoms_for_box(box_id)
      .expect("Box (3, 3, 3) should be present in ForceTaskBoxContainer")
      .collect();

    assert_eq!(particles.len(), 1);

    let particle = particles.first().unwrap();
    let expected_position = Vector3::new(55., 95.0, 55.0);
    assert_eq!(particle.get_position(), expected_position);

    let proxy = particle
      .as_any()
      .downcast_ref::<ParticlePositionProxy>()
      .expect("atoms_for_box should return particle position proxies");
    assert_eq!(
      *proxy.particle_placement(),
      ParticlePlacement {
        x: AxisPlacement::Normal,
        y: AxisPlacement::Normal,
        z: AxisPlacement::Normal,
      }
    );
  }

  #[test]
  fn test_get_needed_box_id_periodic_one_box() {
    let container_config = BoxContainerConfig {
      box_type: FeC,
      box_length: Vector3::new(1., 1., 1.),
      box_count: 100,
      box_count_dim: Vector3::new(10, 10, 10),
      world_size: Vector3::new(10., 10., 10.),
    };

    let box_id = get_id_simulation_box(&Vector3::new(3, 3, 3), &container_config.box_count_dim);

    let mut expected_ids : Vec<usize> = Vec::new();

    for x in 2..=4 {
      for y in 2..=4 {
        for z in 2..=4 {
          let temp_id = get_id_simulation_box(&Vector3::new(x, y, z), &container_config.box_count_dim);
          expected_ids.push(temp_id);
        }
      }
    }

    let obtained_ids = get_needed_box_id_periodic(&vec![box_id], &container_config);

    assert_eq!(obtained_ids.len(), expected_ids.len());

    for id_ in &expected_ids {
      assert!(obtained_ids.contains(id_));
    }
  }

  #[test]
  fn test_get_needed_box_id_periodic_three_boxes() {
    let container_config = BoxContainerConfig {
      box_type: FeC,
      box_length: Vector3::new(1., 1., 1.),
      box_count: 100,
      box_count_dim: Vector3::new(10, 10, 10),
      world_size: Vector3::new(10., 10., 10.),
    };

    let mut argument : Vec<usize> = Vec::new();

    for x in 3..=5 {
      let box_id = get_id_simulation_box(&Vector3::new(x, 3, 3), &container_config.box_count_dim);
      argument.push(box_id);
    }

    let mut temp: HashSet<usize> = HashSet::new();

    for x in 2..=6 {
      for y in 2..=4 {
        for z in 2..=4 {
          let temp_id = get_id_simulation_box(&Vector3::new(x, y, z), &container_config.box_count_dim);
          temp.insert(temp_id);
        }
      }
    }

    let mut expected_ids: Vec<usize> = temp.into_iter().collect();
    expected_ids.sort();

    let mut obtained_ids = get_needed_box_id_periodic(&argument, &container_config);
    obtained_ids.sort();

    assert_eq!(obtained_ids.len(), expected_ids.len());

    for id_ in &expected_ids {
      assert!(obtained_ids.contains(id_));
    }
  }

  #[test]
  fn test_get_needed_box_id_periodic_one_box_left_border() {
    let container_config = BoxContainerConfig {
      box_type: FeC,
      box_length: Vector3::new(1., 1., 1.),
      box_count: 100,
      box_count_dim: Vector3::new(10, 10, 10),
      world_size: Vector3::new(10., 10., 10.),
    };

    let box_id = get_id_simulation_box(&Vector3::new(0, 4, 4), &container_config.box_count_dim);

    let mut expected_ids : Vec<usize> = Vec::new();

    for x in 0..=1 {
      for y in 3..=5 {
        for z in 3..=5 {
          let temp_id = get_id_simulation_box(&Vector3::new(x, y, z), &container_config.box_count_dim);
          expected_ids.push(temp_id);
        }
      }
    }

    for y in 3..=5 {
      for z in 3..=5 {
        let temp_id = get_id_simulation_box(&Vector3::new(9, y, z), &container_config.box_count_dim);
        expected_ids.push(temp_id);
      }
    }

    let obtained_ids = get_needed_box_id_periodic(&vec![box_id], &container_config);

    assert_eq!(obtained_ids.len(), expected_ids.len());

    for id_ in &expected_ids {
      assert!(obtained_ids.contains(id_));
    }
  }
  // FeFe cutoff: R + D = 3.15 + 0.2 = 3.35 = box_size.
  //
  // Box (0,2,2) is the left-x-edge home. Its atoms carry AxisPlacement::Left, so their
  // proxy x = real_x + world_size_x.  Distances that matter (all y,z equal → pure x):
  //
  //   home atom         real x=0.5  → proxy 34.0
  //   box-1 atom        real x=3.5  → proxy (Left) 37.0    distance  3.0  < 3.35  interacts
  //   box-9 atom        real x=33.0 → proxy (Normal) 33.0  distance  1.0  < 3.35  interacts
  //   box-2 atom        real x=7.0  → NOT in neighbour list (2 boxes away)
  //   box-8 atom        real x=30.0 → NOT in neighbour list (2 boxes away)
  //
  // The test verifies:
  //   (a) box-1 and box-9 IDs appear in particles_j  →  one-away interaction
  //   (b) box-2 and box-8 IDs are absent from particles_j  →  two-away exclusion
  //   (c) force on the home atom is nonzero  →  actual physical interaction happens
  #[test]
  fn periodic_x_edge_interacts_with_one_away_not_two_away() {
    use std::collections::HashSet;

    const BOX_SIZE: f64 = 3.35;
    const NX: usize = 10;
    const NY: usize = 5;
    const NZ: usize = 5;
    let world_size = Vector3::new(
      NX as f64 * BOX_SIZE,
      NY as f64 * BOX_SIZE,
      NZ as f64 * BOX_SIZE,
    );
    let box_count_dim = Vector3::new(NX, NY, NZ);

    let container_config = BoxContainerConfig {
      box_type: FeFe,
      box_length: Vector3::new(BOX_SIZE, BOX_SIZE, BOX_SIZE),
      box_count: NX * NY * NZ,
      box_count_dim,
      world_size,
    };

    let edge_condition = EdgeCondition::Periodic {
      split: EdgeCondition::DEFAULT_SPLIT,
      trigger_small_subtask_size: 1,
    };

    let y_center = (2.0 + 0.5) * BOX_SIZE;
    let z_center = (2.0 + 0.5) * BOX_SIZE;

    let home_box_id = get_id_simulation_box(&Vector3::new(0, 2, 2), &box_count_dim);
    let needed_box_ids = get_needed_box_id_periodic(&vec![home_box_id], &container_config);

    // Two-away boxes must not appear among the needed (1-away) boxes.
    let box_2_id = get_id_simulation_box(&Vector3::new(2, 2, 2), &box_count_dim);
    let box_8_id = get_id_simulation_box(&Vector3::new(8, 2, 2), &box_count_dim);
    assert!(!needed_box_ids.contains(&box_2_id), "box-2 must not be in 1-away neighborhood");
    assert!(!needed_box_ids.contains(&box_8_id), "box-8 must not be in 1-away neighborhood");

    let mut atom_factory = SafeAtomFactory::new(1.0, world_size.z);
    let mut container = BoxContainer::new_local(container_config);

    // home: box(0,2,2)
    let home_atom = atom_factory.get_atom(Fe, Vector3::new(0.5, y_center, z_center), Vector3::zeros());
    let home_id = home_atom.get_id();
    container.add_particle(Arc::new(home_atom));

    // 1-away right: box(1,2,2) — proxy x = 3.5+33.5 = 37.0, distance = 3.0
    let atom_1_right = atom_factory.get_atom(Fe, Vector3::new(3.5, y_center, z_center), Vector3::zeros());
    let atom_1_right_id = atom_1_right.get_id();
    container.add_particle(Arc::new(atom_1_right));

    // 1-away left periodic: box(9,2,2) — proxy x = 33.0 (Normal), distance = 1.0
    let atom_1_left = atom_factory.get_atom(Fe, Vector3::new(33.0, y_center, z_center), Vector3::zeros());
    let atom_1_left_id = atom_1_left.get_id();
    container.add_particle(Arc::new(atom_1_left));

    // 2-away right: box(2,2,2) — absent from neighbor list entirely
    let atom_2_right = atom_factory.get_atom(Fe, Vector3::new(7.0, y_center, z_center), Vector3::zeros());
    let atom_2_right_id = atom_2_right.get_id();
    container.add_particle(Arc::new(atom_2_right));

    // 2-away left periodic: box(8,2,2) — absent from neighbor list entirely
    let atom_2_left = atom_factory.get_atom(Fe, Vector3::new(30.0, y_center, z_center), Vector3::zeros());
    let atom_2_left_id = atom_2_left.get_id();
    container.add_particle(Arc::new(atom_2_left));

    let view = container.into_shared().view_select_boxes(&needed_box_ids);
    let force_container = ForceTaskBoxContainer::new(view, edge_condition);

    // Build particles_j the same way handle_force_batch_task does.
    let mut particles_j: Vec<Box<dyn ForceComputationOperations>> =
      force_container.neighbour_atoms_periodic(home_box_id).collect();
    particles_j.extend(force_container.atoms_for_box(home_box_id).unwrap());

    let j_ids: HashSet<usize> = particles_j.iter().map(|p| p.get_id()).collect();

    // (a) one-away atoms must be visible to the home box
    assert!(j_ids.contains(&atom_1_right_id), "box-1 atom must appear in neighbour list");
    assert!(j_ids.contains(&atom_1_left_id), "box-9 (periodic) atom must appear in neighbour list");

    // (b) two-away atoms must be invisible — they never enter the force computation
    assert!(!j_ids.contains(&atom_2_right_id), "box-2 atom must NOT appear in neighbour list");
    assert!(!j_ids.contains(&atom_2_left_id), "box-8 atom must NOT appear in neighbour list");

    // (c) running the actual force computation must produce a nonzero force on the home atom
    let particles_i: Vec<Box<dyn ForceComputationOperations>> =
      force_container.atoms_for_box(home_box_id).unwrap().collect();

    let result = compute_forces_potential(&particles_i, &particles_j);
    let force = result.fp.get(&home_id).expect("home atom must be in force result").force;
    assert!(
      force.magnitude() > 0.0,
      "home atom must have nonzero force from its one-away neighbours (got {:?})",
      force
    );
  }
