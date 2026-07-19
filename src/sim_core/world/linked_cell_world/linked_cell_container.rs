use std::sync::Arc;

use nalgebra::Vector3;

use crate::particle::Particle;
use crate::persistence::dto::world::boxed::box_container::BoxContainerDTO;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::cell::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::{
  SimBoxEdge, get_coordinates_from_simulation_box_id, get_id_simulation_box,
};
use crate::sim_core::world::boxed_world::box_task::force_task_box_container::particle_proxy::{
  AxisPlacement, ParticlePlacement, ParticlePositionProxy, new_particle_position_proxy,
};
use crate::sim_core::world::computation::ForceComputationOperations;
use crate::utils::cube::Cube;

pub struct LinkedCellContainerOld {
  header: Cube<i32>,
  link: Vec<i32>,
  cell: Vec<Vector3<i32>>,
  particles: Vec<Option<Arc<Particle>>>,
  config: BoxContainerConfig,
  edge_condition: EdgeCondition,
}

impl LinkedCellContainerOld {
  pub fn new_empty(num_particles: usize, config: BoxContainerConfig, edge_condition: EdgeCondition) -> Self {
    let (mx, my, mz) = (
      config.box_count_dim.x,
      config.box_count_dim.y,
      config.box_count_dim.z,
    );
    LinkedCellContainerOld {
      header: Cube::new_with_value(mx, my, mz, -1i32),
      link: vec![-1i32; num_particles],
      cell: vec![Vector3::new(-1, -1, -1); num_particles],
      particles: vec![None; num_particles],
      config,
      edge_condition,
    }
  }

  pub fn new(
    particles: Vec<Arc<Particle>>,
    config: BoxContainerConfig,
    edge_condition: EdgeCondition,
  ) -> Self {
    let n = particles.len();
    let (mx, my, mz) = (
      config.box_count_dim.x,
      config.box_count_dim.y,
      config.box_count_dim.z,
    );

    let mut particle_slots: Vec<Option<Arc<Particle>>> = vec![None; n];
    for particle in particles {
      let id = particle.get_id();
      particle_slots[id] = Some(particle);
    }

    LinkedCellContainerOld {
      header: Cube::new_with_value(mx, my, mz, -1i32),
      link: vec![-1i32; n],
      cell: vec![Vector3::new(-1, -1, -1); n],
      particles: particle_slots,
      config,
      edge_condition,
    }
  }

  pub fn add_particle(&mut self, particle: Arc<Particle>) {
    let id = particle.get_id();
    debug_assert!(self.link[id] == -1, "particle {} already sorted", id);
    debug_assert!(self.particles[id].is_none(), "particle {} already present", id);

    let coords = self.config.box_coordinates_for_position(particle.get_position());
    let (kx, ky, kz) = (coords.x, coords.y, coords.z);

    let old_head = *self.header.get(kx, ky, kz).unwrap();
    self.link[id] = old_head;
    self.header.set(kx, ky, kz, id as i32).unwrap();
    self.cell[id] = Vector3::new(coords.x as i32, coords.y as i32, coords.z as i32);
    self.particles[id] = Some(particle);
  }

  /// Assigns each particle to a cell using the linked-cell method.
  ///
  /// After this call:
  /// - `header[kx][ky][kz]` is the index of the last particle placed in that cell (-1 if empty).
  /// - `link[i]` is the index of the previous particle in the same cell as particle i (-1 if none).
  /// - `cell[i]` is the (kx, ky, kz) cell coordinates of particle i.
  pub fn sort(&mut self) {
    self.header.fill(-1);
    for l in self.link.iter_mut() {
      *l = -1;
    }

    for i in 0..self.particles.len() {
      if let Some(particle) = &self.particles[i] {
        let coords = self.config.box_coordinates_for_position(particle.get_position());
        let (kx, ky, kz) = (coords.x, coords.y, coords.z);

        let old_head = *self.header.get(kx, ky, kz).unwrap();
        self.link[i] = old_head;
        self.header.set(kx, ky, kz, i as i32).unwrap();
        self.cell[i] = Vector3::new(coords.x as i32, coords.y as i32, coords.z as i32);
      }
    }
  }

  /// Returns an iterator over all particles in the given cell (identified by flat cell id).
  pub fn particles_in_cell(&self, cell_id: usize) -> impl Iterator<Item = Arc<Particle>> {
    let coords =
      get_coordinates_from_simulation_box_id(cell_id, &self.config.box_count_dim);
    let mut idx = *self.header.get(coords.x, coords.y, coords.z).unwrap_or(&-1);

    let mut result = Vec::with_capacity(4);
    while idx >= 0 {
      let particle = self.particles[idx as usize].as_ref().unwrap();
      result.push(Arc::clone(particle));
      idx = self.link[idx as usize];
    }
    result.into_iter()
  }

  /// Returns an iterator over particles in the given cell wrapped in a position proxy
  /// that accounts for periodic boundary conditions.
  pub fn atoms_for_cell(
    &self,
    cell_id: usize,
  ) -> impl Iterator<Item = ParticlePositionProxy> {
    let coords =
      get_coordinates_from_simulation_box_id(cell_id, &self.config.box_count_dim);
    let cell_count = self.config.box_count_dim;
    let world_size = self.config.world_size;

    let x_edge = if coords.x == 0 {
      SimBoxEdge::LeftEdge
    } else if coords.x == cell_count.x - 1 {
      SimBoxEdge::RightEdge
    } else {
      SimBoxEdge::Normal
    };

    let y_edge = if coords.y == 0 {
      SimBoxEdge::LeftEdge
    } else if coords.y == cell_count.y - 1 {
      SimBoxEdge::RightEdge
    } else {
      SimBoxEdge::Normal
    };

    let z_edge = if coords.z == 0 {
      SimBoxEdge::LeftEdge
    } else if coords.z == cell_count.z - 1 {
      SimBoxEdge::RightEdge
    } else {
      SimBoxEdge::Normal
    };

    let x_axis_placement = match (self.edge_condition, x_edge) {
      (EdgeCondition::Simple { .. }, _) => unimplemented!("not yet"),
      (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge) => {
        AxisPlacement::Left
      }
      (
        EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll,
        SimBoxEdge::Normal | SimBoxEdge::RightEdge,
      ) => AxisPlacement::Normal,
    };

    let y_axis_placement = match (self.edge_condition, y_edge) {
      (EdgeCondition::Simple { .. }, _) => unimplemented!("not yet"),
      (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge) => {
        AxisPlacement::Left
      }
      (
        EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll,
        SimBoxEdge::Normal | SimBoxEdge::RightEdge,
      ) => AxisPlacement::Normal,
    };

    let z_axis_placement = match (self.edge_condition, z_edge) {
      (EdgeCondition::Simple { .. }, _) => unimplemented!("not yet"),
      (EdgeCondition::Periodic { .. }, _) => AxisPlacement::Normal,
      (EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge) => AxisPlacement::Left,
      (EdgeCondition::PeriodicAll, SimBoxEdge::Normal | SimBoxEdge::RightEdge) => {
        AxisPlacement::Normal
      }
    };

    let placement = ParticlePlacement {
      x: x_axis_placement,
      y: y_axis_placement,
      z: z_axis_placement,
    };

    self
      .particles_in_cell(cell_id)
      .map(move |particle| new_particle_position_proxy(particle, placement, world_size))
  }

  /// Returns an iterator over all particles in the 26 cells surrounding the given cell,
  /// wrapped in position proxies that account for periodic boundary conditions.
  pub fn neighbour_atoms_periodic(
    &self,
    cell_id: usize,
  ) -> impl Iterator<Item = ParticlePositionProxy> {
    let coords = get_coordinates_from_simulation_box_id(cell_id, &self.config.box_count_dim);
    let cell_count = self.config.box_count_dim;
    let world_size = self.config.world_size;

    let x_edge = if coords.x == 0 {
      SimBoxEdge::LeftEdge
    } else if coords.x == cell_count.x - 1 {
      SimBoxEdge::RightEdge
    } else {
      SimBoxEdge::Normal
    };
    let y_edge = if coords.y == 0 {
      SimBoxEdge::LeftEdge
    } else if coords.y == cell_count.y - 1 {
      SimBoxEdge::RightEdge
    } else {
      SimBoxEdge::Normal
    };
    let z_edge = if coords.z == 0 {
      SimBoxEdge::LeftEdge
    } else if coords.z == cell_count.z - 1 {
      SimBoxEdge::RightEdge
    } else {
      SimBoxEdge::Normal
    };

    let x = coords.x as isize;
    let y = coords.y as isize;
    let z = coords.z as isize;

    let mut result = Vec::new();

    for x_offset in -1..=1isize {
      for y_offset in -1..=1isize {
        for z_offset in -1..=1isize {
          if x_offset == 0 && y_offset == 0 && z_offset == 0 {
            continue;
          }

          let x_axis_placement = match (self.edge_condition, x_edge, x_offset) {
            (EdgeCondition::Simple { .. }, _, _) => AxisPlacement::Normal,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge, 0 | 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::RightEdge, 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, _, _) => AxisPlacement::Normal,
          };

          let y_axis_placement = match (self.edge_condition, y_edge, y_offset) {
            (EdgeCondition::Simple { .. }, _, _) => AxisPlacement::Normal,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge, 0 | 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, SimBoxEdge::RightEdge, 1) => AxisPlacement::Left,
            (EdgeCondition::Periodic { .. } | EdgeCondition::PeriodicAll, _, _) => AxisPlacement::Normal,
          };

          let z_axis_placement = match (self.edge_condition, z_edge, z_offset) {
            (EdgeCondition::Simple { .. }, _, _) => AxisPlacement::Normal,
            (EdgeCondition::Periodic { .. }, _, _) => AxisPlacement::Normal,
            (EdgeCondition::PeriodicAll, SimBoxEdge::LeftEdge, 0 | 1) => AxisPlacement::Left,
            (EdgeCondition::PeriodicAll, SimBoxEdge::RightEdge, 1) => AxisPlacement::Left,
            (EdgeCondition::PeriodicAll, _, _) => AxisPlacement::Normal,
          };

          let nx = (((x + x_offset) + cell_count.x as isize) as usize) % cell_count.x;
          let ny = (((y + y_offset) + cell_count.y as isize) as usize) % cell_count.y;
          let nz = (((z + z_offset) + cell_count.z as isize) as usize) % cell_count.z;

          let neighbour_id = get_id_simulation_box(&Vector3::new(nx, ny, nz), &cell_count);
          let placement = ParticlePlacement {
            x: x_axis_placement,
            y: y_axis_placement,
            z: z_axis_placement,
          };

          for particle in self.particles_in_cell(neighbour_id) {
            result.push(new_particle_position_proxy(particle, placement, world_size));
          }
        }
      }
    }

    result.into_iter()
  }

  pub fn header(&self) -> &Cube<i32> {
    &self.header
  }

  pub fn link(&self) -> &[i32] {
    &self.link
  }

  pub fn cell(&self) -> &[Vector3<i32>] {
    &self.cell
  }

  pub fn particles(&self) -> &[Option<Arc<Particle>>] {
    &self.particles
  }

  pub fn particle_mut(&mut self, id: usize) -> &mut Particle {
    Arc::make_mut(self.particles[id].as_mut().unwrap())
  }

  pub fn config(&self) -> &BoxContainerConfig {
    &self.config
  }

  pub fn edge_condition(&self) -> EdgeCondition {
    self.edge_condition
  }

  pub fn to_transfer_struct(&self) -> BoxContainerDTO {
    let atoms = self.particles.iter()
      .filter_map(|p| p.as_ref())
      .map(|p| p.to_transfer_struct())
      .collect();
    BoxContainerDTO { atoms, config: self.config }
  }
}

#[cfg(test)]
mod tests {
  use super::*;
  use crate::data::constants::ATOMIC_MASS_FE;
  use crate::data::types::{AtomType, InteractionType};
  use crate::particle::atom::Atom;
  use crate::persistence::json::particle_config::read_particle_config_from_json_str;
  use crate::sim_core::world::cell::box_container_config::new_config;
  use crate::sim_core::world::boxed_world::box_container::sim_box::get_id_simulation_box;
  use crate::sim_core::world::boxed_world::box_task::force_task_box_container::particle_proxy::ParticlePositionProxy;
  use std::collections::HashSet;

  fn make_config() -> BoxContainerConfig {
    // 10×10×10 world, 2×2×2 grid → each cell is 5×5×5
    BoxContainerConfig {
      box_type: InteractionType::FeFe,
      box_length: Vector3::new(5.0, 5.0, 5.0),
      box_count: 8,
      box_count_dim: Vector3::new(2, 2, 2),
      world_size: Vector3::new(10.0, 10.0, 10.0),
    }
  }

  fn make_4x4x4_config() -> BoxContainerConfig {
    // 20×20×20 world, 4×4×4 grid → each cell is 5×5×5
    // Cells at indices 1 and 2 (per axis) are interior (not at any edge).
    BoxContainerConfig {
      box_type: InteractionType::FeFe,
      box_length: Vector3::new(5.0, 5.0, 5.0),
      box_count: 64,
      box_count_dim: Vector3::new(4, 4, 4),
      world_size: Vector3::new(20.0, 20.0, 20.0),
    }
  }

  fn make_container_one_per_cell(
    config: BoxContainerConfig,
    edge_condition: EdgeCondition,
  ) -> LinkedCellContainerOld {
    let cell_count = config.box_count_dim;
    let cell_size = config.box_length;
    let mut particles = Vec::new();
    let mut id = 0;
    for x in 0..cell_count.x {
      for y in 0..cell_count.y {
        for z in 0..cell_count.z {
          let pos = Vector3::new(
            (x as f64 + 0.5) * cell_size.x,
            (y as f64 + 0.5) * cell_size.y,
            (z as f64 + 0.5) * cell_size.z,
          );
          particles.push(make_particle(id, pos));
          id += 1;
        }
      }
    }
    let mut container = LinkedCellContainerOld::new(particles, config, edge_condition);
    container.sort();
    container
  }

  fn make_particle(id: usize, position: Vector3<f64>) -> Arc<Particle> {
    Arc::new(Particle::Atom(Atom::new(
      id,
      AtomType::Fe,
      ATOMIC_MASS_FE,
      position,
      Vector3::zeros(),
      Vector3::zeros(),
      Vector3::zeros(),
      0.0,
    )))
  }

  #[test]
  fn sort_single_particle() {
    // Particle at (1, 1, 1) → cell (0, 0, 0)
    let p0 = make_particle(0, Vector3::new(1.0, 1.0, 1.0));
    let mut container =
      LinkedCellContainerOld::new(vec![p0], make_config(), EdgeCondition::PeriodicAll);

    container.sort();

    assert_eq!(*container.header().get(0, 0, 0).unwrap(), 0);
    assert_eq!(container.link()[0], -1);
    assert_eq!(container.cell()[0], Vector3::new(0, 0, 0));
    // All other cells remain empty
    assert_eq!(*container.header().get(1, 1, 1).unwrap(), -1);
  }

  #[test]
  fn sort_two_particles_different_cells() {
    // p0 at (1, 1, 1) → cell (0, 0, 0)
    // p1 at (6, 6, 6) → cell (1, 1, 1)
    let p0 = make_particle(0, Vector3::new(1.0, 1.0, 1.0));
    let p1 = make_particle(1, Vector3::new(6.0, 6.0, 6.0));
    let mut container =
      LinkedCellContainerOld::new(vec![p0, p1], make_config(), EdgeCondition::PeriodicAll);

    container.sort();

    assert_eq!(*container.header().get(0, 0, 0).unwrap(), 0);
    assert_eq!(*container.header().get(1, 1, 1).unwrap(), 1);
    assert_eq!(container.link()[0], -1);
    assert_eq!(container.link()[1], -1);
    assert_eq!(container.cell()[0], Vector3::new(0, 0, 0));
    assert_eq!(container.cell()[1], Vector3::new(1, 1, 1));
  }

  #[test]
  fn sort_multiple_particles_mixed_cells() {
    // p0 (1,1,1) → cell (0,0,0)
    // p1 (6,1,1) → cell (1,0,0)
    // p2 (2,1,1) → cell (0,0,0)  — chains onto p0
    // p3 (7,7,7) → cell (1,1,1)
    // p4 (3,2,2) → cell (0,0,0)  — chains onto p2
    let particles = vec![
      make_particle(0, Vector3::new(1.0, 1.0, 1.0)),
      make_particle(1, Vector3::new(6.0, 1.0, 1.0)),
      make_particle(2, Vector3::new(2.0, 1.0, 1.0)),
      make_particle(3, Vector3::new(7.0, 7.0, 7.0)),
      make_particle(4, Vector3::new(3.0, 2.0, 2.0)),
    ];
    let mut container =
      LinkedCellContainerOld::new(particles, make_config(), EdgeCondition::PeriodicAll);

    container.sort();

    // cell (0,0,0): p4 → p2 → p0 → -1
    assert_eq!(*container.header().get(0, 0, 0).unwrap(), 4);
    assert_eq!(container.link()[4], 2);
    assert_eq!(container.link()[2], 0);
    assert_eq!(container.link()[0], -1);

    // cell (1,0,0): only p1
    assert_eq!(*container.header().get(1, 0, 0).unwrap(), 1);
    assert_eq!(container.link()[1], -1);

    // cell (1,1,1): only p3
    assert_eq!(*container.header().get(1, 1, 1).unwrap(), 3);
    assert_eq!(container.link()[3], -1);

    // cell (0,1,0) untouched
    assert_eq!(*container.header().get(0, 1, 0).unwrap(), -1);

    // cell assignments
    assert_eq!(container.cell()[0], Vector3::new(0, 0, 0));
    assert_eq!(container.cell()[1], Vector3::new(1, 0, 0));
    assert_eq!(container.cell()[2], Vector3::new(0, 0, 0));
    assert_eq!(container.cell()[3], Vector3::new(1, 1, 1));
    assert_eq!(container.cell()[4], Vector3::new(0, 0, 0));
  }

  #[test]
  fn sort_two_particles_same_cell() {
    // p0 at (1, 1, 1) → cell (0, 0, 0)
    // p1 at (2, 2, 2) → cell (0, 0, 0)
    // p1 is added last, so header points to p1; link[1] = 0, link[0] = -1
    let p0 = make_particle(0, Vector3::new(1.0, 1.0, 1.0));
    let p1 = make_particle(1, Vector3::new(2.0, 2.0, 2.0));
    let mut container =
      LinkedCellContainerOld::new(vec![p0, p1], make_config(), EdgeCondition::PeriodicAll);

    container.sort();

    assert_eq!(*container.header().get(0, 0, 0).unwrap(), 1);
    assert_eq!(container.link()[1], 0);
    assert_eq!(container.link()[0], -1);
    assert_eq!(container.cell()[0], Vector3::new(0, 0, 0));
    assert_eq!(container.cell()[1], Vector3::new(0, 0, 0));
  }

  #[test]
  fn particles_in_cell_empty() {
    // All particles land in cell (0,0,0); querying cell (1,0,0) must yield nothing.
    let particles = vec![
      make_particle(0, Vector3::new(1.0, 1.0, 1.0)),
      make_particle(1, Vector3::new(2.0, 2.0, 2.0)),
    ];
    let mut container =
      LinkedCellContainerOld::new(particles, make_config(), EdgeCondition::PeriodicAll);
    container.sort();

    let cell_id = get_id_simulation_box(&Vector3::new(1, 0, 0), &make_config().box_count_dim);
    let result: Vec<_> = container.particles_in_cell(cell_id).collect();

    assert!(result.is_empty());
  }

  #[test]
  fn particles_in_cell_one_particle() {
    // p0 in cell (0,0,0), p1 in cell (1,0,0) — query cell (1,0,0) returns only p1.
    let p0 = make_particle(0, Vector3::new(1.0, 1.0, 1.0));
    let p1 = make_particle(1, Vector3::new(6.0, 1.0, 1.0));
    let mut container =
      LinkedCellContainerOld::new(vec![p0, p1], make_config(), EdgeCondition::PeriodicAll);
    container.sort();

    let cell_id = get_id_simulation_box(&Vector3::new(1, 0, 0), &make_config().box_count_dim);
    let result: Vec<_> = container.particles_in_cell(cell_id).collect();

    assert_eq!(result.len(), 1);
    assert_eq!(result[0].get_id(), 1);
  }

  #[test]
  fn particles_in_cell_many_particles() {
    // p0, p1, p2 all in cell (0,0,0); p3 in cell (1,1,1).
    // Linked list inserts at head, so the iterator yields p2 → p1 → p0.
    let particles = vec![
      make_particle(0, Vector3::new(1.0, 1.0, 1.0)),
      make_particle(1, Vector3::new(2.0, 2.0, 2.0)),
      make_particle(2, Vector3::new(3.0, 3.0, 3.0)),
      make_particle(3, Vector3::new(7.0, 7.0, 7.0)),
    ];
    let mut container =
      LinkedCellContainerOld::new(particles, make_config(), EdgeCondition::PeriodicAll);
    container.sort();

    let cell_id = get_id_simulation_box(&Vector3::new(0, 0, 0), &make_config().box_count_dim);
    let result: Vec<_> = container.particles_in_cell(cell_id).collect();

    assert_eq!(result.len(), 3);
    assert_eq!(result[0].get_id(), 2);
    assert_eq!(result[1].get_id(), 1);
    assert_eq!(result[2].get_id(), 0);
  }

  #[test]
  fn sort_four_particles_same_cell() {
    // All four particles fall into cell (0,0,0).
    // They are added in index order, so the chain runs: header=3 → 2 → 1 → 0 → -1
    let particles = vec![
      make_particle(0, Vector3::new(1.0, 1.0, 1.0)),
      make_particle(1, Vector3::new(2.0, 2.0, 2.0)),
      make_particle(2, Vector3::new(3.0, 3.0, 3.0)),
      make_particle(3, Vector3::new(4.0, 4.0, 4.0)),
    ];
    let mut container =
      LinkedCellContainerOld::new(particles, make_config(), EdgeCondition::PeriodicAll);

    container.sort();

    assert_eq!(*container.header().get(0, 0, 0).unwrap(), 3);
    assert_eq!(container.link()[3], 2);
    assert_eq!(container.link()[2], 1);
    assert_eq!(container.link()[1], 0);
    assert_eq!(container.link()[0], -1);

    for i in 0..4 {
      assert_eq!(container.cell()[i], Vector3::new(0, 0, 0));
    }
    // all other cells are empty
    assert_eq!(*container.header().get(1, 0, 0).unwrap(), -1);
    assert_eq!(*container.header().get(0, 1, 0).unwrap(), -1);
    assert_eq!(*container.header().get(1, 1, 1).unwrap(), -1);
  }

  // ── atoms_for_cell ──────────────────────────────────────────────────────────

  #[test]
  fn atoms_for_cell_interior_normal_placement() {
    // Cell (2,2,2) in a 4×4×4 grid is interior — not at any edge.
    // The proxy must leave the position unchanged and report all-Normal placement.
    let config = make_4x4x4_config();
    let cell_id = get_id_simulation_box(&Vector3::new(2, 2, 2), &config.box_count_dim);
    let container = make_container_one_per_cell(
      config,
      EdgeCondition::Periodic {
        trigger_small_subtask_size: 1,
        split: EdgeCondition::DEFAULT_SPLIT,
      },
    );

    let particles: Vec<_> = container.atoms_for_cell(cell_id).collect();

    assert_eq!(particles.len(), 1);
    assert_eq!(particles[0].get_position(), Vector3::new(12.5, 12.5, 12.5));
    let proxy = particles[0]
      .as_any()
      .downcast_ref::<ParticlePositionProxy>()
      .expect("atoms_for_cell should return ParticlePositionProxy");
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
  fn atoms_for_cell_left_x_edge_left_placement() {
    // Cell (0,2,2) sits at the left x-edge.
    // AxisPlacement::Left shifts the proxy x by +world_size.x (20): 2.5 → 22.5.
    let config = make_4x4x4_config();
    let cell_id = get_id_simulation_box(&Vector3::new(0, 2, 2), &config.box_count_dim);
    let container = make_container_one_per_cell(
      config,
      EdgeCondition::Periodic {
        trigger_small_subtask_size: 1,
        split: EdgeCondition::DEFAULT_SPLIT,
      },
    );

    let particles: Vec<_> = container.atoms_for_cell(cell_id).collect();

    assert_eq!(particles.len(), 1);
    assert_eq!(particles[0].get_position(), Vector3::new(22.5, 12.5, 12.5));
    let proxy = particles[0]
      .as_any()
      .downcast_ref::<ParticlePositionProxy>()
      .expect("atoms_for_cell should return ParticlePositionProxy");
    assert_eq!(
      *proxy.particle_placement(),
      ParticlePlacement {
        x: AxisPlacement::Left,
        y: AxisPlacement::Normal,
        z: AxisPlacement::Normal,
      }
    );
  }

  // ── neighbour_atoms_periodic ─────────────────────────────────────────────────

  #[test]
  fn neighbour_atoms_periodic_interior_all_normal() {
    // Center cell (2,2,2) in a 4×4×4 grid.
    // None of the 3×3×3 neighbourhood crosses a periodic boundary from the center
    // cell's perspective, so every neighbour gets all-Normal placement.
    let config = make_4x4x4_config();
    let cell_id = get_id_simulation_box(&Vector3::new(2, 2, 2), &config.box_count_dim);
    let container = make_container_one_per_cell(
      config,
      EdgeCondition::Periodic {
        trigger_small_subtask_size: 1,
        split: EdgeCondition::DEFAULT_SPLIT,
      },
    );

    let particles: Vec<_> = container.neighbour_atoms_periodic(cell_id).collect();

    assert_eq!(particles.len(), 26);
    for particle in &particles {
      let proxy = particle
        .as_any()
        .downcast_ref::<ParticlePositionProxy>()
        .expect("neighbour_atoms_periodic should return ParticlePositionProxy");
      assert_eq!(
        *proxy.particle_placement(),
        ParticlePlacement {
          x: AxisPlacement::Normal,
          y: AxisPlacement::Normal,
          z: AxisPlacement::Normal,
        }
      );
    }
  }

  #[test]
  fn neighbour_atoms_periodic_left_x_edge_mixed_placement() {
    // Center cell (0,2,2): x_edge = LeftEdge under EdgeCondition::Periodic.
    // x_offset=-1 wraps to x=3 → Normal x-placement (9 particles)
    // x_offset= 0 stays at x=0, excluding center → Left x-placement (8 particles)
    // x_offset=+1 moves to x=1 → Left x-placement (9 particles)
    // y and z are interior (index 2 of 4), so y/z placement is always Normal.
    let config = make_4x4x4_config();
    let cell_id = get_id_simulation_box(&Vector3::new(0, 2, 2), &config.box_count_dim);
    let container = make_container_one_per_cell(
      config,
      EdgeCondition::Periodic {
        trigger_small_subtask_size: 1,
        split: EdgeCondition::DEFAULT_SPLIT,
      },
    );

    let particles: Vec<_> = container.neighbour_atoms_periodic(cell_id).collect();

    assert_eq!(particles.len(), 26);

    let mut normal_x_count = 0usize;
    let mut left_x_count = 0usize;
    for particle in &particles {
      let proxy = particle
        .as_any()
        .downcast_ref::<ParticlePositionProxy>()
        .expect("neighbour_atoms_periodic should return ParticlePositionProxy");
      assert_eq!(proxy.particle_placement().y, AxisPlacement::Normal);
      assert_eq!(proxy.particle_placement().z, AxisPlacement::Normal);
      match proxy.particle_placement().x {
        AxisPlacement::Normal => normal_x_count += 1,
        AxisPlacement::Left => left_x_count += 1,
      }
    }

    assert_eq!(normal_x_count, 9);  // neighbours at x=3 (wrapped by x_offset=-1)
    assert_eq!(left_x_count, 17);   // neighbours at x=0 (8) and x=1 (9)
  }

  fn make_10x10x10_config() -> BoxContainerConfig {
    // 100×100×100 world, 10×10×10 grid → each cell is 10×10×10
    BoxContainerConfig {
      box_type: InteractionType::FeFe,
      box_length: Vector3::new(10.0, 10.0, 10.0),
      box_count: 1000,
      box_count_dim: Vector3::new(10, 10, 10),
      world_size: Vector3::new(100.0, 100.0, 100.0),
    }
  }

  #[test]
  fn neighbour_atoms_periodic_left_x_edge_explicit_positions_and_placements() {
    // Mirrors test_force_task_box_container_with_one_particle_in_each_box_center_neighbour_atoms_periodic_left_edge_placement.
    // Grid 10×10×10, cells 10×10×10, world 100×100×100. One particle per cell at cell center.
    // Center cell (0,3,3): x=LeftEdge, y=Normal, z=Normal under EdgeCondition::Periodic.
    //
    // x_offset=-1 → wraps to x=9, raw pos_x=95  → Normal x (9 particles)
    // x_offset= 0 → stays  at x=0, raw pos_x=5  → Left x, shifted to 105 (8 particles, center excluded)
    // x_offset=+1 → moves  to x=1, raw pos_x=15 → Left x, shifted to 115 (9 particles)
    let config = make_10x10x10_config();
    let cell_id = get_id_simulation_box(&Vector3::new(0, 3, 3), &config.box_count_dim);
    let container = make_container_one_per_cell(
      config,
      EdgeCondition::Periodic {
        trigger_small_subtask_size: 1,
        split: EdgeCondition::DEFAULT_SPLIT,
      },
    );

    let particles: Vec<_> = container.neighbour_atoms_periodic(cell_id).collect();

    assert_eq!(particles.len(), 26);

    let normal = ParticlePlacement { x: AxisPlacement::Normal, y: AxisPlacement::Normal, z: AxisPlacement::Normal };
    let left_x = ParticlePlacement { x: AxisPlacement::Left,   y: AxisPlacement::Normal, z: AxisPlacement::Normal };

    let mut expected: Vec<(Vector3<f64>, ParticlePlacement)> = vec![
      // x_offset=-1 → x=9, pos_x=95, Normal (all 9 y/z combos)
      (Vector3::new(95.0, 25.0, 25.0), normal),
      (Vector3::new(95.0, 25.0, 35.0), normal),
      (Vector3::new(95.0, 25.0, 45.0), normal),
      (Vector3::new(95.0, 35.0, 25.0), normal),
      (Vector3::new(95.0, 35.0, 35.0), normal),
      (Vector3::new(95.0, 35.0, 45.0), normal),
      (Vector3::new(95.0, 45.0, 25.0), normal),
      (Vector3::new(95.0, 45.0, 35.0), normal),
      (Vector3::new(95.0, 45.0, 45.0), normal),
      // x_offset=0 → x=0, pos_x=5+100=105, Left (center excluded → 8)
      (Vector3::new(105.0, 25.0, 25.0), left_x),
      (Vector3::new(105.0, 25.0, 35.0), left_x),
      (Vector3::new(105.0, 25.0, 45.0), left_x),
      (Vector3::new(105.0, 35.0, 25.0), left_x),
      (Vector3::new(105.0, 35.0, 45.0), left_x),
      (Vector3::new(105.0, 45.0, 25.0), left_x),
      (Vector3::new(105.0, 45.0, 35.0), left_x),
      (Vector3::new(105.0, 45.0, 45.0), left_x),
      // x_offset=+1 → x=1, pos_x=15+100=115, Left (all 9)
      (Vector3::new(115.0, 25.0, 25.0), left_x),
      (Vector3::new(115.0, 25.0, 35.0), left_x),
      (Vector3::new(115.0, 25.0, 45.0), left_x),
      (Vector3::new(115.0, 35.0, 25.0), left_x),
      (Vector3::new(115.0, 35.0, 35.0), left_x),
      (Vector3::new(115.0, 35.0, 45.0), left_x),
      (Vector3::new(115.0, 45.0, 25.0), left_x),
      (Vector3::new(115.0, 45.0, 35.0), left_x),
      (Vector3::new(115.0, 45.0, 45.0), left_x),
    ];

    for particle in &particles {
      let proxy = particle
        .as_any()
        .downcast_ref::<ParticlePositionProxy>()
        .expect("neighbour_atoms_periodic should return ParticlePositionProxy");
      let pos = particle.get_position();
      let placement = *proxy.particle_placement();
      let idx = expected
        .iter()
        .position(|(ep, epl)| *ep == pos && *epl == placement)
        .expect("particle has unexpected position or placement");
      expected.remove(idx);
    }

    assert!(expected.is_empty());
  }

  // ── insertion-order equivalence (fixture-based) ─────────────────────────────

  /// Verifies that building a LinkedCellContainerOld via `new()` + `sort()` (with
  /// particles pre-sorted by id) produces the same cell contents as building it
  /// via `new_empty()` + `add_particle()` (with particles in the shuffled order
  /// found in the JSON fixture).  The two construction paths use different array
  /// indexing strategies, so this test checks whether insertion order matters.
  #[test]
  fn insertion_order_does_not_affect_cells() {
    let fixture_path = concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/particles_initial.json");
    let fixture = std::fs::read_to_string(fixture_path).expect("fixture file not found");

    let particle_config = read_particle_config_from_json_str(&fixture)
      .expect("failed to parse fixture");

    // Positions are now in unitless (Angstrom) coordinates: ~0.8..25.3.
    // A 30×30×30 world comfortably contains all particles.
    let particles = particle_config.atoms;
    let world_size = Vector3::new(30.0, 30.0, 30.0);
    let config = new_config(&particles, world_size);
    let edge_condition = EdgeCondition::PeriodicAll;

    let num_particles = particles.iter().map(|p| p.get_id()).max().unwrap() + 1;

    // Container A: particles sorted by id → new() + sort()
    let mut sorted = particles.clone();
    let arcs_sorted: Vec<Arc<Particle>> = sorted.iter().map(|p| Arc::new(p.clone())).collect();
    let mut container_a = LinkedCellContainerOld::new(arcs_sorted, config, edge_condition);
    container_a.sort();

    // Container B: original JSON order (shuffled ids) → new_empty() + add_particle()
    let mut container_b = LinkedCellContainerOld::new_empty(num_particles, config, edge_condition);
    for particle in &particles {
      container_b.add_particle(Arc::new(particle.clone()));
    }

    let pos_bits = |v: Vector3<f64>| -> [u64; 3] {
      [v.x.to_bits(), v.y.to_bits(), v.z.to_bits()]
    };

    for cell_id in 0..config.box_count {
      let ids_a: HashSet<usize> = container_a.particles_in_cell(cell_id)
        .map(|p| p.get_id())
        .collect();
      let ids_b: HashSet<usize> = container_b.particles_in_cell(cell_id)
        .map(|p| p.get_id())
        .collect();
      assert_eq!(ids_a, ids_b, "cell {cell_id}: particles_in_cell id sets differ");

      let atoms_a: HashSet<(usize, [u64; 3])> = container_a.atoms_for_cell(cell_id)
        .map(|p| (p.get_id(), pos_bits(p.get_position())))
        .collect();
      let atoms_b: HashSet<(usize, [u64; 3])> = container_b.atoms_for_cell(cell_id)
        .map(|p| (p.get_id(), pos_bits(p.get_position())))
        .collect();
      assert_eq!(atoms_a, atoms_b, "cell {cell_id}: atoms_for_cell (id, pos) sets differ");

      let nbrs_a: HashSet<(usize, [u64; 3])> = container_a.neighbour_atoms_periodic(cell_id)
        .map(|p| (p.get_id(), pos_bits(p.get_position())))
        .collect();
      let nbrs_b: HashSet<(usize, [u64; 3])> = container_b.neighbour_atoms_periodic(cell_id)
        .map(|p| (p.get_id(), pos_bits(p.get_position())))
        .collect();
      assert_eq!(nbrs_a, nbrs_b, "cell {cell_id}: neighbour_atoms_periodic (id, pos) sets differ");
    }
  }

  // ── ID-index invariant (particles[i].get_id() == i) ─────────────────────────

  /// `new()` stores particles in insertion order, not by particle ID.
  /// With the shuffled fixture the invariant `particles[i].id == i` does NOT hold.
  /// This test documents that bug: it will fail until `new()` is fixed to
  /// index by ID, or callers are required to pre-sort.
  #[test]
  fn new_with_unsorted_particles_id_matches_index() {
    let fixture_path =
      concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/particles_initial.json");
    let fixture = std::fs::read_to_string(fixture_path).expect("fixture file not found");
    let particle_config =
      read_particle_config_from_json_str(&fixture).expect("failed to parse fixture");

    let particles = particle_config.atoms;
    let world_size = Vector3::new(30.0, 30.0, 30.0);
    let config = new_config(&particles, world_size);

    let arcs: Vec<Arc<Particle>> = particles.iter().map(|p| Arc::new(p.clone())).collect();
    let mut container = LinkedCellContainerOld::new(arcs, config, EdgeCondition::PeriodicAll);
    container.sort();

    for (i, slot) in container.particles().iter().enumerate() {
      if let Some(particle) = slot {
        assert_eq!(
          particle.get_id(),
          i,
          "particles[{i}] holds particle with id {}",
          particle.get_id()
        );
      }
    }
  }

  /// `new_empty()` + `add_particle()` stores each particle at index equal to its ID,
  /// so the invariant `particles[i].id == i` holds regardless of insertion order.
  #[test]
  fn new_empty_add_particle_id_matches_index() {
    let fixture_path =
      concat!(env!("CARGO_MANIFEST_DIR"), "/tests/fixtures/particles_initial.json");
    let fixture = std::fs::read_to_string(fixture_path).expect("fixture file not found");
    let particle_config =
      read_particle_config_from_json_str(&fixture).expect("failed to parse fixture");

    let particles = particle_config.atoms;
    let world_size = Vector3::new(30.0, 30.0, 30.0);
    let config = new_config(&particles, world_size);

    let num_particles = particles.iter().map(|p| p.get_id()).max().unwrap() + 1;
    let mut container =
      LinkedCellContainerOld::new_empty(num_particles, config, EdgeCondition::PeriodicAll);
    for particle in &particles {
      container.add_particle(Arc::new(particle.clone()));
    }

    for (i, slot) in container.particles().iter().enumerate() {
      if let Some(particle) = slot {
        assert_eq!(
          particle.get_id(),
          i,
          "particles[{i}] holds particle with id {}",
          particle.get_id()
        );
      }
    }
  }
}
