use std::collections::HashSet;
use std::sync::Arc;

use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::{
  SimBoxEdge, SimulationBox, get_coordinates_from_simulation_box_id, get_id_simulation_box,
};
use crate::sim_core::world::boxed_world::box_task::force_task_box_container::particle_proxy::{
  AxisPlacement, ParticlePlacement, ParticlePositionProxy, new_particle_position_proxy,
};
use crate::sim_core::world::computation::ForceComputationOperations;
use nalgebra::Vector3;

pub mod particle_proxy;

pub fn get_needed_box_id_periodic(box_ids: &Vec<usize>, config: &BoxContainerConfig) -> Vec<usize> {
  let mut temp: HashSet<usize> = HashSet::new();
  let box_count_dim = config.box_count_dim;

  for id in box_ids {
    let coordinates = get_coordinates_from_simulation_box_id(*id, &box_count_dim);

    for x_offset in -1..=1isize {
      for y_offset in -1..=1isize {
        for z_offset in -1..=1isize {
          let new_x = ((coordinates.x as isize + x_offset + box_count_dim.x as isize) as usize)
            % box_count_dim.x;
          let new_y = ((coordinates.y as isize + y_offset + box_count_dim.y as isize) as usize)
            % box_count_dim.y;
          let new_z = ((coordinates.z as isize + z_offset + box_count_dim.z as isize) as usize)
            % box_count_dim.z;

          let new_id = get_id_simulation_box(&Vector3::new(new_x, new_y, new_z), &box_count_dim);
          temp.insert(new_id);
        }
      }
    }
  }

  temp.into_iter().collect()
}

pub struct ForceTaskBoxContainer {
  container: BoxContainer<Option<Arc<SimulationBox>>>,
  edge_condition: EdgeCondition,
}

impl ForceTaskBoxContainer {
  pub fn new(
    container: BoxContainer<Option<Arc<SimulationBox>>>,
    edge_condition: EdgeCondition,
  ) -> Self {
    ForceTaskBoxContainer {
      container,
      edge_condition,
    }
  }

  pub fn atoms_for_box(
    &self,
    box_id: usize,
  ) -> Option<Box<dyn Iterator<Item = Box<dyn ForceComputationOperations>>>> {
    let config = self.container.config();
    let world_size = config.world_size;

    let sim_box = self.container.get_box(box_id)?;

    let x_axis_placement = match (self.edge_condition, sim_box.sim_box_placement().x) {
      (EdgeCondition::Simple { .. }, _) => unimplemented!("not yet"),
      (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge) => AxisPlacement::Left,
      (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::Normal | SimBoxEdge::RightEdge) => {
        AxisPlacement::Normal
      }
    };

    let y_axis_placement = match (self.edge_condition, sim_box.sim_box_placement().y) {
      (EdgeCondition::Simple { .. }, _) => unimplemented!("not yet"),
      (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge) => AxisPlacement::Left,
      (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::Normal | SimBoxEdge::RightEdge) => {
        AxisPlacement::Normal
      }
    };

    let z_axis_placement = match (self.edge_condition, sim_box.sim_box_placement().z) {
      (EdgeCondition::Simple { .. }, _) => unimplemented!("not yet"),
      (EdgeCondition::Periodic { .. }, _) => AxisPlacement::Normal,
      (EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge) => AxisPlacement::Left,
      (EdgeCondition::PeriodicAll, SimBoxEdge::Normal | SimBoxEdge::RightEdge) => AxisPlacement::Normal,
    };

    let placement = ParticlePlacement {
      x: x_axis_placement,
      y: y_axis_placement,
      z: z_axis_placement
    };

    let proxies: Vec<Box<dyn ForceComputationOperations>> = sim_box
      .particles()
      .values()
      .map(|particle| {
        Box::new(new_particle_position_proxy(
          Arc::clone(particle),
          placement,
          world_size,
        )) as Box<dyn ForceComputationOperations>
      })
      .collect();

    Some(Box::new(proxies.into_iter()))
  }

  pub fn neighbour_atoms_periodic(
    &self,
    box_id: usize,
  ) -> Box<dyn Iterator<Item = Box<dyn ForceComputationOperations>>> {
    if matches!(self.edge_condition, EdgeCondition::Simple { .. }) {
      panic!("this method doesn't work correctly for simple edge condition.")
    }

    let config = self.container.config();
    let world_size = config.world_size;
    let box_count_dim = config.box_count_dim;

    let coordinates = get_coordinates_from_simulation_box_id(box_id, &box_count_dim);
    let sim_box_placement = self
      .container
      .get_box(box_id)
      .expect("Target box must be present in ForceTaskBoxContainer")
      .sim_box_placement();

    let x = coordinates.x as isize;
    let y = coordinates.y as isize;
    let z = coordinates.z as isize;

    let mut proxies: Vec<Box<dyn ForceComputationOperations>> = Vec::new();

    for x_offset in -1..=1isize {
      for y_offset in -1..=1isize {
        for z_offset in -1..=1isize {
          if x_offset == 0 && y_offset == 0 && z_offset == 0 {
            continue;
          }

          let x_axis_placement = match (self.edge_condition, sim_box_placement.x, x_offset) {
            (EdgeCondition::Simple { .. }, _, _) => AxisPlacement::Normal,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge, 0 | 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::RightEdge, 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, _, _) => AxisPlacement::Normal,
          };

          let y_axis_placement = match (self.edge_condition, sim_box_placement.y, y_offset) {
            (EdgeCondition::Simple { .. }, _, _) => AxisPlacement::Normal,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge, 0 | 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::RightEdge, 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, _, _) => AxisPlacement::Normal,
          };

          let z_axis_placement = match (self.edge_condition, sim_box_placement.z, z_offset) {
            (EdgeCondition::Simple { .. }, _, _) => AxisPlacement::Normal,
            (EdgeCondition::Periodic { .. }, _, _) => AxisPlacement::Normal,
            (EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge, 0 | 1) => AxisPlacement::Left,
            (EdgeCondition::PeriodicAll, SimBoxEdge::RightEdge, 1) => AxisPlacement::Left,
            (EdgeCondition::PeriodicAll, _, _) => AxisPlacement::Normal,
          };

          let x_ = (((x + x_offset) + box_count_dim.x as isize) as usize) % box_count_dim.x;
          let y_ = (((y + y_offset) + box_count_dim.y as isize) as usize) % box_count_dim.y;
          let z_ = (((z + z_offset) + box_count_dim.z as isize) as usize) % box_count_dim.z;

          let neighbour_id = get_id_simulation_box(&Vector3::new(x_, y_, z_), &box_count_dim);

          let sim_box = self
            .container
            .get_box(neighbour_id)
            .expect("Neighbour box must be present in ForceTaskBoxContainer");
          let placement = ParticlePlacement {
            x: x_axis_placement,
            y: y_axis_placement,
            z: z_axis_placement
          };
          for particle in sim_box.particles().values() {
            proxies.push(Box::new(new_particle_position_proxy(
              Arc::clone(particle),
              placement,
              world_size,
            )));
          }
        }
      }
    }

    Box::new(proxies.into_iter())
  }

  pub fn view(&self) -> &BoxContainer<Option<Arc<SimulationBox>>> {
    &self.container
  }
}

#[cfg(test)]
#[path = "force_task_box_container_test.rs"]
mod tests;