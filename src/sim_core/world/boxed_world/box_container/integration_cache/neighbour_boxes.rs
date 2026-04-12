use std::sync::Arc;
use nalgebra::Vector3;
use crate::data::types::AtomType;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::{get_coordinates_from_simulation_box_id, SimBoxEdge};
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::ForceComputationOperations;

#[derive(Copy, Clone)]
pub enum AxisPlacement {
  Left,
  Normal
}

#[derive(Copy, Clone)]
pub struct ParticlePlacement {
  x: AxisPlacement,
  y: AxisPlacement,
}

pub struct ParticlePositionProxy {
  particle: Arc<Particle>,
  particle_placement: ParticlePlacement,
  world_size: Vector3<f64>,
}

pub fn new_particle_position_proxy(particle: Arc<Particle>, particle_placement: ParticlePlacement,
                                   world_size: Vector3<f64>) -> ParticlePositionProxy {
  ParticlePositionProxy {
    particle,
    particle_placement,
    world_size,
  }
}

impl ForceComputationOperations for ParticlePositionProxy {
  fn get_id(&self) -> usize {
    self.particle.get_id() as usize
  }

  fn get_position(&self) -> Vector3<f64> {
    let mut position = self.particle.get_position().clone();

    position.x = match self.particle_placement.x {
      AxisPlacement::Left => position.x + self.world_size.x,
      AxisPlacement::Normal => position.x,
    };

    position.y = match self.particle_placement.y {
      AxisPlacement::Left => position.y + self.world_size.y,
      AxisPlacement::Normal => position.y,
    };

    position
  }

  fn get_type(&self) -> AtomType {
    self.particle.get_type()
  }

  fn get_mass(&self) -> f64 {
    self.particle.get_mass()
  }

  fn prototype_clone(&self) -> Box<dyn ForceComputationOperations> {
    self.particle.prototype_clone()
  }
}

impl BoxContainer {
  pub fn atoms_for_force_computation_of_single_integration_box(&self, box_id: usize)
                                                               -> Box<dyn Iterator<Item = Box<dyn ForceComputationOperations>> + '_> {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    let world_size = self.container_size;
    let simulation_box = self.integration_box_cache.get(coordinates.x, coordinates.y, coordinates.z).unwrap();

    Box::new(simulation_box.particles().values().map(move |particle| {
      let x_axis_placement = match (self.edge_condition(), simulation_box.sim_box_placement().x) {
        (EdgeCondition::Simple, _) => AxisPlacement::Normal,
        (EdgeCondition::Periodic, SimBoxEdge::LeftEdge) => AxisPlacement::Left,
        (EdgeCondition::Periodic, SimBoxEdge::Normal | SimBoxEdge::RightEdge) => AxisPlacement::Normal,
      };

      let y_axis_placement = match (self.edge_condition(), simulation_box.sim_box_placement().y) {
        (EdgeCondition::Simple, _) => AxisPlacement::Normal,
        (EdgeCondition::Periodic, SimBoxEdge::LeftEdge) => AxisPlacement::Left,
        (EdgeCondition::Periodic, SimBoxEdge::Normal | SimBoxEdge::RightEdge) => AxisPlacement::Normal,
      };

      let proxy = ParticlePositionProxy {
        particle: Arc::clone(&particle),
        particle_placement: ParticlePlacement {
          x: x_axis_placement,
          y: y_axis_placement,
        },
        world_size,
      };

      Box::new(proxy) as Box<dyn ForceComputationOperations>
    }))
  }

  pub fn integration_particles_of_neighbour_boxes_simple(&self, box_id: usize) -> impl Iterator<Item = Arc<Particle>> {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    let box_count_dim = self.box_count_dim();

    let x = coordinates.x as isize;
    let y = coordinates.y as isize;
    let z = coordinates.z as isize;

    (-1..=1).flat_map(move |x_offset| {
      (-1..=1).flat_map(move |y_offset| {
        (-1..=1).filter_map(move |z_offset| {
          if x_offset == 0 && y_offset == 0 && z_offset == 0 {
            return None;
          }
          let x_ = x + x_offset;
          let y_ = y + y_offset;
          let z_ = z + z_offset;

          if x_ >= 0 && x_ < box_count_dim.x as isize &&
            y_ >= 0 && y_ < box_count_dim.y as isize &&
            z_ >= 0 && z_ < box_count_dim.z as isize {
            self.integration_box_cache.get(x_ as usize, y_ as usize, z_ as usize)
              .map(|sim_box| sim_box.particles().values().cloned())
          } else {
            None
          }
        }).flatten()
      })
    })
  }

  pub fn integration_particles_of_neighbour_boxes_periodic(&self, box_id: usize) -> Box<dyn Iterator<Item = Box<dyn ForceComputationOperations>> + '_> {
    let world_size = self.container_size().clone();
    let box_count_dim = self.box_count_dim();

    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    let simulation_box_placement = self.integration_box_cache
      .get(coordinates.x, coordinates.y, coordinates.z).unwrap().sim_box_placement();

    let x = coordinates.x as isize;
    let y = coordinates.y as isize;
    let z = coordinates.z as isize;

    Box::new((-1..=1).flat_map(move |x_offset| {
      (-1..=1).flat_map(move |y_offset| {
        (-1..=1).filter_map(move |z_offset| {
          if x_offset == 0 && y_offset == 0 && z_offset == 0 {
            return None;
          }

          let x_axis_placement = match (self.edge_condition, simulation_box_placement.x, x_offset) {
            (EdgeCondition::Simple, _ ,_) => AxisPlacement::Normal,
            (EdgeCondition::Periodic, SimBoxEdge::LeftEdge, 0 | 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic, SimBoxEdge::RightEdge, 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic, _, _) => AxisPlacement::Normal,
          };

          let y_axis_placement = match (self.edge_condition, simulation_box_placement.y, y_offset) {
            (EdgeCondition::Simple, _ ,_) => AxisPlacement::Normal,
            (EdgeCondition::Periodic, SimBoxEdge::LeftEdge, 0 | 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic, SimBoxEdge::RightEdge, 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic, _, _) => AxisPlacement::Normal,
          };

          let x_ = (((x + x_offset) + box_count_dim.x as isize) as usize) % box_count_dim.x;
          let y_ = (((y + y_offset) + box_count_dim.y as isize) as usize) % box_count_dim.y;
          let z_ = (((z + z_offset) + box_count_dim.z as isize) as usize) % box_count_dim.z;

          self.integration_box_cache.get(x_, y_, z_)
            .map(|sim_box| sim_box.particles().values()
              .map(move |particle| {
              let proxy = ParticlePositionProxy {
                particle: Arc::clone(&particle),
                particle_placement: ParticlePlacement {
                  x: x_axis_placement,
                  y: y_axis_placement,
                },
                world_size,
              };
              Box::new(proxy) as Box<dyn ForceComputationOperations>
              })
            )
        }).flatten()
      })
    }))
  }

  pub fn integration_particles_of_neighbour_boxes(&self, box_id: usize) -> Box<dyn Iterator<Item = Box<dyn ForceComputationOperations>> + '_> {
    let world_size = self.container_size;
    match self.edge_condition() {
      EdgeCondition::Simple => {
        Box::new(
          self.integration_particles_of_neighbour_boxes_simple(box_id)
            .map(move |particle| {
              Box::new(ParticlePositionProxy {
                particle: Arc::clone(&particle),
                particle_placement: ParticlePlacement {
                  x: AxisPlacement::Normal,
                  y: AxisPlacement::Normal,
                },
                world_size,
              }) as Box<dyn ForceComputationOperations>
            })
        )
      },
      EdgeCondition::Periodic => {
        self.integration_particles_of_neighbour_boxes_periodic(box_id)
      }
    }
  }
}