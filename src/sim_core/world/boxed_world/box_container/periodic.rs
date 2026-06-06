use std::collections::HashSet;

use nalgebra::Vector3;

use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::{
  get_coordinates_from_simulation_box_id, get_id_simulation_box,
};

/// Periodic offset of box coordinates; wraps into `[0, dim)`.
pub fn periodic_box_coordinates(
  coordinates: Vector3<usize>,
  offset: Vector3<isize>,
  box_count_dim: &Vector3<usize>,
) -> Vector3<usize> {
  let nx = ((coordinates.x as isize + offset.x + box_count_dim.x as isize) as usize) % box_count_dim.x;
  let ny = ((coordinates.y as isize + offset.y + box_count_dim.y as isize) as usize) % box_count_dim.y;
  let nz = ((coordinates.z as isize + offset.z + box_count_dim.z as isize) as usize) % box_count_dim.z;
  Vector3::new(nx, ny, nz)
}

/// 3×3×3 periodic halo around `box_ids` (includes each center box).
pub fn get_needed_box_id_periodic(box_ids: &[usize], config: &BoxContainerConfig) -> Vec<usize> {
  let mut temp: HashSet<usize> = HashSet::new();
  let box_count_dim = config.box_count_dim;

  for &id in box_ids {
    let coordinates = get_coordinates_from_simulation_box_id(id, &box_count_dim);

    for x_offset in -1..=1isize {
      for y_offset in -1..=1isize {
        for z_offset in -1..=1isize {
          let new_coords = periodic_box_coordinates(
            coordinates,
            Vector3::new(x_offset, y_offset, z_offset),
            &box_count_dim,
          );
          let new_id = get_id_simulation_box(&new_coords, &box_count_dim);
          temp.insert(new_id);
        }
      }
    }
  }

  temp.into_iter().collect()
}

/// Face-adjacent periodic neighbours only (6-connected), excluding `box_id`.
pub fn face_neighbor_box_ids_periodic(box_id: usize, box_count_dim: &Vector3<usize>) -> Vec<usize> {
  let coordinates = get_coordinates_from_simulation_box_id(box_id, box_count_dim);
  let offsets = [
    Vector3::new(-1, 0, 0),
    Vector3::new(1, 0, 0),
    Vector3::new(0, -1, 0),
    Vector3::new(0, 1, 0),
    Vector3::new(0, 0, -1),
    Vector3::new(0, 0, 1),
  ];

  offsets
    .into_iter()
    .map(|offset| {
      let coords = periodic_box_coordinates(coordinates, offset, box_count_dim);
      get_id_simulation_box(&coords, box_count_dim)
    })
    .collect()
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::data::types::InteractionType::FeC;

  fn sample_config() -> BoxContainerConfig {
    BoxContainerConfig {
      box_type: FeC,
      box_length: Vector3::new(1., 1., 1.),
      box_count: 100,
      box_count_dim: Vector3::new(10, 10, 10),
      world_size: Vector3::new(10., 10., 10.),
    }
  }

  #[test]
  fn get_needed_box_id_periodic_one_box_has_27_cells() {
    let container_config = sample_config();
    let box_id = get_id_simulation_box(&Vector3::new(3, 3, 3), &container_config.box_count_dim);

    let mut expected_ids: Vec<usize> = Vec::new();
    for x in 2..=4 {
      for y in 2..=4 {
        for z in 2..=4 {
          expected_ids.push(get_id_simulation_box(
            &Vector3::new(x, y, z),
            &container_config.box_count_dim,
          ));
        }
      }
    }

    let mut obtained_ids = get_needed_box_id_periodic(&[box_id], &container_config);
    obtained_ids.sort();
    expected_ids.sort();

    assert_eq!(obtained_ids.len(), 27);
    assert_eq!(obtained_ids, expected_ids);
  }

  #[test]
  fn face_neighbor_box_ids_periodic_returns_six() {
    let container_config = sample_config();
    let box_id = get_id_simulation_box(&Vector3::new(3, 3, 3), &container_config.box_count_dim);
    let neighbors = face_neighbor_box_ids_periodic(box_id, &container_config.box_count_dim);
    assert_eq!(neighbors.len(), 6);
    assert!(!neighbors.contains(&box_id));
  }
}
