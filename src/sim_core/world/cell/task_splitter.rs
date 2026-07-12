use std::collections::HashMap;
use std::sync::Arc;

use nalgebra::Vector3;

use crate::sim_core::world::boxed_world::box_container::sim_box::get_id_simulation_box;
use crate::sim_core::world::cell::box_container_config::BoxContainerConfig;

#[derive(Debug, Clone, Copy)]
pub enum TaskSplitVariant {
  Floor,
  FloorBox { x: usize, y: usize },
}

pub struct TaskSplitter {
  variant: TaskSplitVariant,
}

impl TaskSplitter {
  pub fn new(variant: TaskSplitVariant) -> Self {
    TaskSplitter { variant }
  }

  pub fn split(
    &self,
    num_of_tasks: usize,
    config: &BoxContainerConfig,
  ) -> HashMap<usize, Arc<Vec<usize>>> {
    match self.variant {
      TaskSplitVariant::Floor => Self::split_floor(num_of_tasks, config),
      TaskSplitVariant::FloorBox { x, y } => Self::split_boxes(num_of_tasks, config, x, y),
    }
  }

  pub fn split_floor(
    num_of_tasks: usize,
    config: &BoxContainerConfig,
  ) -> HashMap<usize, Arc<Vec<usize>>> {
    let nx = config.box_count_dim.x;
    let ny = config.box_count_dim.y;
    let nz = config.box_count_dim.z;

    let z_per_task = nz / num_of_tasks;
    let remainder = nz % num_of_tasks;

    let mut mapping = HashMap::with_capacity(num_of_tasks);
    let mut z_start = 0;

    for task_id in 0..num_of_tasks {
      let z_count = z_per_task + if task_id < remainder { 1 } else { 0 };
      let z_end = z_start + z_count;

      let cell_ids = (z_start..z_end)
        .flat_map(|z| {
          (0..ny).flat_map(move |y| {
            (0..nx).map(move |x| {
              get_id_simulation_box(&Vector3::new(x, y, z), &config.box_count_dim)
            })
          })
        })
        .collect();

      mapping.insert(task_id, Arc::new(cell_ids));
      z_start = z_end;
    }

    mapping
  }

  pub fn split_boxes(
    num_of_tasks: usize,
    config: &BoxContainerConfig,
    x_splits: usize,
    y_splits: usize,
  ) -> HashMap<usize, Arc<Vec<usize>>> {
    let nx = config.box_count_dim.x;
    let ny = config.box_count_dim.y;
    let nz = config.box_count_dim.z;

    let num_of_tasks_per_floor = x_splits * y_splits;
    let num_of_floors = num_of_tasks / num_of_tasks_per_floor;

    let x_per_box = nx / x_splits;
    let x_remainder = nx % x_splits;

    let y_per_box = ny / y_splits;
    let y_remainder = ny % y_splits;

    let z_per_floor = nz / num_of_floors;
    let z_remainder = nz % num_of_floors;

    let mut mapping = HashMap::with_capacity(num_of_floors * num_of_tasks_per_floor);
    let mut task_id = 0;
    let mut z_start = 0;

    for floor_id in 0..num_of_floors {
      let z_count = z_per_floor + if floor_id < z_remainder { 1 } else { 0 };
      let z_end = z_start + z_count;

      let mut y_start = 0;
      for y_box_id in 0..y_splits {
        let y_count = y_per_box + if y_box_id < y_remainder { 1 } else { 0 };
        let y_end = y_start + y_count;

        let mut x_start = 0;
        for x_box_id in 0..x_splits {
          let x_count = x_per_box + if x_box_id < x_remainder { 1 } else { 0 };
          let x_end = x_start + x_count;

          let cell_ids = (z_start..z_end)
            .flat_map(|z| {
              (y_start..y_end).flat_map(move |y| {
                (x_start..x_end).map(move |x| {
                  get_id_simulation_box(&Vector3::new(x, y, z), &config.box_count_dim)
                })
              })
            })
            .collect();

          mapping.insert(task_id, Arc::new(cell_ids));
          task_id += 1;
          x_start = x_end;
        }
        y_start = y_end;
      }
      z_start = z_end;
    }

    mapping
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::data::types::InteractionType;
  use std::ops::Range;

  fn make_config(box_count_dim: Vector3<usize>) -> BoxContainerConfig {
    BoxContainerConfig {
      box_type: InteractionType::FeFe,
      box_length: Vector3::new(1.0, 1.0, 1.0),
      box_count: box_count_dim.x * box_count_dim.y * box_count_dim.z,
      box_count_dim,
      world_size: Vector3::new(
        box_count_dim.x as f64,
        box_count_dim.y as f64,
        box_count_dim.z as f64,
      ),
    }
  }

  fn expected_cell_ids(
    x_range: Range<usize>,
    y_range: Range<usize>,
    z_range: Range<usize>,
    box_count_dim: &Vector3<usize>,
  ) -> Vec<usize> {
    let mut ids = Vec::new();
    for z in z_range.clone() {
      for y in y_range.clone() {
        for x in x_range.clone() {
          ids.push(get_id_simulation_box(&Vector3::new(x, y, z), box_count_dim));
        }
      }
    }
    ids
  }

  // 7x4x5 grid (nx, ny, nz), split with FloorBox { x: 3, y: 2 }, num_of_tasks = 15.
  //
  // num_of_tasks_per_floor = 3 * 2 = 6, num_of_floors = 15 / 6 = 2 (integer division,
  // so only 12 of the requested 15 tasks are actually produced).
  //
  // x: nx=7 / 3 splits -> per_box=2 remainder=1 -> ranges [0,3) [3,5) [5,7)
  // y: ny=4 / 2 splits -> per_box=2 remainder=0 -> ranges [0,2) [2,4)
  // z: nz=5 / 2 floors -> per_floor=2 remainder=1 -> ranges [0,3) [3,5)
  #[test]
  fn split_boxes_matches_hand_computed_ranges_for_7x4x5_grid() {
    let box_count_dim = Vector3::new(7, 4, 5);
    let config = make_config(box_count_dim);

    let mapping = TaskSplitter::split_boxes(15, &config, 3, 2);

    assert_eq!(mapping.len(), 12);

    let expected_ranges: Vec<(usize, Range<usize>, Range<usize>, Range<usize>)> = vec![
      (0, 0..3, 0..2, 0..3),
      (1, 3..5, 0..2, 0..3),
      (2, 5..7, 0..2, 0..3),
      (3, 0..3, 2..4, 0..3),
      (4, 3..5, 2..4, 0..3),
      (5, 5..7, 2..4, 0..3),
      (6, 0..3, 0..2, 3..5),
      (7, 3..5, 0..2, 3..5),
      (8, 5..7, 0..2, 3..5),
      (9, 0..3, 2..4, 3..5),
      (10, 3..5, 2..4, 3..5),
      (11, 5..7, 2..4, 3..5),
    ];

    for (task_id, x_range, y_range, z_range) in expected_ranges {
      let expected = expected_cell_ids(x_range, y_range, z_range, &box_count_dim);
      let actual = mapping
        .get(&task_id)
        .unwrap_or_else(|| panic!("missing task_id {}", task_id));
      assert_eq!(actual.as_ref(), &expected, "task_id {} cell ids mismatch", task_id);
    }

    let mut all_ids: Vec<usize> = mapping.values().flat_map(|v| v.iter().copied()).collect();
    all_ids.sort_unstable();
    let expected_all: Vec<usize> = (0..(7 * 4 * 5)).collect();
    assert_eq!(all_ids, expected_all, "every cell must be assigned exactly once");
  }
}
