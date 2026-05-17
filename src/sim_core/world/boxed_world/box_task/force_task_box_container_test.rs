use crate::data::InteractionType::FeC;
  use crate::data::types::AtomType::C;
  use crate::particle::SafeAtomFactory;
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
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic);

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
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic);

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
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic);
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
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic);
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
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic);

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
    let force_container = ForceTaskBoxContainer::new(view, EdgeCondition::Periodic);

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