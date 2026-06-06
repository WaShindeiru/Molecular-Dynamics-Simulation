use std::collections::HashSet;
use std::sync::Arc;

use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::{
  SimulationBox, get_coordinates_from_simulation_box_id, get_id_simulation_box,
};
use nalgebra::Vector3;
use crate::utils::cube::Cube;

impl BoxContainer<Option<Arc<SimulationBox>>> {
  pub fn get_box(&self, box_id: usize) -> Option<Arc<SimulationBox>> {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
    self
      .simulation_boxes
      .get(coordinates.x, coordinates.y, coordinates.z)
      .unwrap()
      .clone()
  }

  /// Ensures each `box_id` is present as `Some(empty box)` in the force view.
  pub fn ensure_boxes_exist(&mut self, box_ids: &[usize]) {
    let template = BoxContainer::<SimulationBox>::new_local(self.config);
    for &box_id in box_ids {
      if self.get_box(box_id).is_some() {
        continue;
      }
      let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
      let empty_box = Arc::new(template.get_box(box_id).clone());
      *self
        .simulation_boxes
        .get_mut(coordinates.x, coordinates.y, coordinates.z)
        .unwrap() = Some(empty_box);
    }
  }

  /// Merges working primaries with a static halo force view.
  /// Returns `None` when the two sources share any particle id.
  /// TODO: this can be optimized
  pub fn merge_force_views(
    mut working: BoxContainer<SimulationBox>,
    static_halo: &Self,
  ) -> Option<Self> {
    assert_eq!(working.config().box_count_dim, static_halo.config().box_count_dim);

    let working_ids: HashSet<usize> = working.box_id_cache().keys().copied().collect();
    let static_ids = static_halo_particle_ids(static_halo);
    if working_ids.iter().any(|id| static_ids.contains(id)) {
      return None;
    }

    for ((x, y, z), static_cell) in static_halo.simulation_boxes.iter_with_coords() {
      let Some(static_box) = static_cell.as_ref() else {
        continue;
      };
      let box_id =
        get_id_simulation_box(&Vector3::new(x, y, z), &static_halo.config.box_count_dim);
      for particle in static_box.particles().values() {
        working.add_particle_with_box_id(Arc::clone(particle), box_id);
      }
    }

    Some(to_force_option_view(working))
  }
}

fn static_halo_particle_ids(
  static_halo: &BoxContainer<Option<Arc<SimulationBox>>>,
) -> HashSet<usize> {
  static_halo
    .simulation_boxes
    .iter()
    .filter_map(|cell| cell.as_ref())
    .flat_map(|sim_box| sim_box.particles().keys().copied())
    .collect()
}

fn to_force_option_view(
  working: BoxContainer<SimulationBox>,
) -> BoxContainer<Option<Arc<SimulationBox>>> {
  let config = working.config;
  let box_id_cache = working.box_id_cache;
  let (nx, ny, nz) = working.simulation_boxes.dimensions();
  let mut result: Cube<Option<Arc<SimulationBox>>> = Cube::new(nx, ny, nz);

  for ((x, y, z), sim_box) in working.simulation_boxes.iter_with_coords() {
    let value = if sim_box.empty() {
      None
    } else {
      Some(Arc::new(sim_box.clone()))
    };
    result.set(x, y, z, value).unwrap();
  }

  BoxContainer {
    config,
    simulation_boxes: result,
    box_id_cache,
  }
}

#[cfg(test)]
mod merge_force_views_test;
