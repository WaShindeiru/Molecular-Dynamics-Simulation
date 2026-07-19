use nalgebra::Vector3;

use crate::data::types::AtomType;
use crate::particle::Particle;
use crate::particle::atom::Atom;
use crate::persistence::dto::world::boxed::box_container::BoxContainerDTO;
use crate::sim_core::world::boundary_constraint::EdgeCondition;
use crate::sim_core::world::boxed_world::box_container::sim_box::{
  get_coordinates_from_simulation_box_id, get_id_simulation_box,
};
use crate::sim_core::world::cell::box_container_config::BoxContainerConfig;
use crate::sim_core::world::cell::fixed_position_particle::FixedPositionParticle;
use crate::utils::cube::Cube;

/// Never read before being overwritten: every slot's `present` flag starts `false`, and
/// nothing in this module reads a slot without checking `present` (`is_complete`) first.
fn placeholder_particle() -> Particle {
  Particle::Atom(Atom::new(
    0,
    AtomType::Fe,
    0.,
    Vector3::zeros(),
    Vector3::zeros(),
    Vector3::zeros(),
    Vector3::zeros(),
    0.,
  ))
}

/// Resolves one axis of a neighbour-cell offset: wraps it (periodic) or bounds-checks it
/// (rigid), returning the neighbour's coordinate on this axis and the shift (-1, 0, or +1,
/// to be multiplied by that axis's world size) that a particle's position needs so pairwise
/// distances to the home cell stay correct across the wraparound. `None` means this offset
/// falls outside the box and (for a non-wrapping axis) the neighbour cell doesn't exist.
fn wrapped_axis(coord: usize, offset: isize, count: usize, wrap: bool) -> Option<(usize, f64)> {
  let raw = coord as isize + offset;

  if wrap {
    let wrapped = ((raw % count as isize) + count as isize) % count as isize;
    let shift = if raw < 0 {
      -1.0
    } else if raw >= count as isize {
      1.0
    } else {
      0.0
    };
    Some((wrapped as usize, shift))
  } else if raw < 0 || raw >= count as isize {
    None
  } else {
    Some((raw as usize, 0.0))
  }
}

/// Per-axis wrap flags for a given edge condition: `Simple` never wraps, `Periodic` wraps
/// x/y only (z stays rigid), `PeriodicAll` wraps all three.
fn wrap_flags(edge_condition: EdgeCondition) -> (bool, bool, bool) {
  match edge_condition {
    EdgeCondition::Simple { .. } => (false, false, false),
    EdgeCondition::Periodic { .. } => (true, true, false),
    EdgeCondition::PeriodicAll => (true, true, true),
  }
}

/// Linked-cell particle container, always backed by a dense `Vec<Particle>` (no per-element
/// `Arc`, no `Option`). While a container is being assembled from asynchronous task results
/// (`new_empty` + `add_particle`), unfilled slots hold a throwaway placeholder and `present`
/// tracks which ids have real data; `is_complete` is the single completeness check callers
/// use before trusting the container for reads (force computation, sorting into cells, etc.).
#[derive(Clone)]
pub struct LinkedCellContainer {
  header: Cube<i32>,
  link: Vec<i32>,
  cell: Vec<Vector3<i32>>,
  particles: Vec<Particle>,
  present: Vec<bool>,
  config: BoxContainerConfig,
  edge_condition: EdgeCondition,
}

impl LinkedCellContainer {
  pub fn new_empty(num_particles: usize, config: BoxContainerConfig, edge_condition: EdgeCondition) -> Self {
    let (mx, my, mz) = (
      config.box_count_dim.x,
      config.box_count_dim.y,
      config.box_count_dim.z,
    );
    LinkedCellContainer {
      header: Cube::new_with_value(mx, my, mz, -1i32),
      link: vec![-1i32; num_particles],
      cell: vec![Vector3::new(-1, -1, -1); num_particles],
      particles: vec![placeholder_particle(); num_particles],
      present: vec![false; num_particles],
      config,
      edge_condition,
    }
  }

  /// Builds a container directly from a complete particle list (does not place particles into
  /// cells — call `sort()` afterward if cell membership is needed).
  pub fn new(particles: Vec<Particle>, config: BoxContainerConfig, edge_condition: EdgeCondition) -> Self {
    let n = particles.len();
    let (mx, my, mz) = (
      config.box_count_dim.x,
      config.box_count_dim.y,
      config.box_count_dim.z,
    );

    let mut slots = vec![placeholder_particle(); n];
    let mut present = vec![false; n];
    for particle in particles {
      let id = particle.get_id();
      slots[id] = particle;
      present[id] = true;
    }

    LinkedCellContainer {
      header: Cube::new_with_value(mx, my, mz, -1i32),
      link: vec![-1i32; n],
      cell: vec![Vector3::new(-1, -1, -1); n],
      particles: slots,
      present,
      config,
      edge_condition,
    }
  }

  pub fn add_particle(&mut self, particle: Particle) {
    let id = particle.get_id();
    debug_assert!(self.link[id] == -1, "particle {} already sorted", id);
    debug_assert!(!self.present[id], "particle {} already present", id);

    let coords = self.config.box_coordinates_for_position(particle.get_position());
    let (kx, ky, kz) = (coords.x, coords.y, coords.z);

    let old_head = *self.header.get(kx, ky, kz).unwrap();
    self.link[id] = old_head;
    self.header.set(kx, ky, kz, id as i32).unwrap();
    self.cell[id] = Vector3::new(coords.x as i32, coords.y as i32, coords.z as i32);
    self.particles[id] = particle;
    self.present[id] = true;
  }

  /// Whether every id has received real particle data. The single check that replaces
  /// converting to a separate dense representation: once this is `true`, every read method
  /// below can be trusted directly, no per-access presence checks needed.
  pub fn is_complete(&self) -> bool {
    self.present.iter().all(|&p| p)
  }

  /// Produces a blank container that keeps only the per-particle metadata that never changes
  /// across iterations (id, type, mass, ...) via `Particle::reset_clone`, with `header`/`link`/
  /// `cell`/`present` at their empty defaults and `config`/`edge_condition` preserved. Meant to
  /// be computed once per simulation run and then reused (via `Clone`) as the starting point for
  /// each iteration's scratch containers, since it doesn't depend on anything that changes
  /// iteration to iteration (position, velocity, force, ... are stripped by `reset_clone`).
  pub fn reset_clone(&self) -> LinkedCellContainer {
    let (mx, my, mz) = (
      self.config.box_count_dim.x,
      self.config.box_count_dim.y,
      self.config.box_count_dim.z,
    );
    let n = self.particles.len();
    LinkedCellContainer {
      header: Cube::new_with_value(mx, my, mz, -1i32),
      link: vec![-1i32; n],
      cell: vec![Vector3::new(-1, -1, -1); n],
      particles: self.particles.iter().map(|p| p.reset_clone()).collect(),
      present: vec![false; n],
      config: self.config,
      edge_condition: self.edge_condition,
    }
  }

  /// Moves the particle at `id` to `new_position` and links it into its new cell. Unlike
  /// `add_particle`, the slot at `id` is expected to already hold a particle (e.g. from
  /// `reset_clone`) rather than a placeholder — only position/cell bookkeeping changes.
  pub fn change_position(&mut self, id: usize, new_position: Vector3<f64>) {
    debug_assert!(self.link[id] == -1, "particle {} already sorted", id);
    debug_assert!(!self.present[id], "particle {} already present", id);

    self.particles[id].update_position(new_position);

    let coords = self.config.box_coordinates_for_position(&new_position);
    let (kx, ky, kz) = (coords.x, coords.y, coords.z);

    let old_head = *self.header.get(kx, ky, kz).unwrap();
    self.link[id] = old_head;
    self.header.set(kx, ky, kz, id as i32).unwrap();
    self.cell[id] = Vector3::new(coords.x as i32, coords.y as i32, coords.z as i32);
    self.present[id] = true;
  }

  /// Assigns each present particle to a cell using the linked-cell method.
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
      if self.present[i] {
        let coords = self.config.box_coordinates_for_position(self.particles[i].get_position());
        let (kx, ky, kz) = (coords.x, coords.y, coords.z);

        let old_head = *self.header.get(kx, ky, kz).unwrap();
        self.link[i] = old_head;
        self.header.set(kx, ky, kz, i as i32).unwrap();
        self.cell[i] = Vector3::new(coords.x as i32, coords.y as i32, coords.z as i32);
      }
    }
  }

  /// Returns the ids of all particles in the given cell (identified by flat cell id).
  pub fn ids_in_cell(&self, cell_id: usize) -> Vec<usize> {
    let coords = get_coordinates_from_simulation_box_id(cell_id, &self.config.box_count_dim);
    let mut idx = *self.header.get(coords.x, coords.y, coords.z).unwrap_or(&-1);

    let mut result = Vec::with_capacity(4);
    while idx >= 0 {
      result.push(idx as usize);
      idx = self.link[idx as usize];
    }
    result
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

  pub fn particles(&self) -> &[Particle] {
    &self.particles
  }

  pub fn particle_mut(&mut self, id: usize) -> &mut Particle {
    &mut self.particles[id]
  }

  pub fn config(&self) -> &BoxContainerConfig {
    &self.config
  }

  pub fn edge_condition(&self) -> EdgeCondition {
    self.edge_condition
  }

  pub fn to_transfer_struct(&self) -> BoxContainerDTO {
    let atoms = self
      .particles
      .iter()
      .zip(self.present.iter())
      .filter_map(|(p, &present)| present.then(|| p.to_transfer_struct()))
      .collect();
    BoxContainerDTO { atoms, config: self.config }
  }

  /// Returns `{id, position}` for every particle in the given cell. Positions are never
  /// shifted here (a cell's own particles are always at zero offset from themselves).
  pub fn atoms_for_cell_fixed(&self, cell_id: usize) -> Vec<FixedPositionParticle> {
    self
      .ids_in_cell(cell_id)
      .into_iter()
      .map(|id| FixedPositionParticle { id, position: *self.particles[id].get_position() })
      .collect()
  }

  /// Returns `{id, position}` for every particle in the (up to 26) neighbour cells around
  /// the given cell, with positions shifted for wrapped axes so pairwise distances to the
  /// home cell are correct. Which axes wrap depends on `edge_condition`: `Simple` wraps
  /// none (out-of-bounds neighbour cells are skipped), `Periodic` wraps x/y only (z is
  /// skipped like `Simple` when out of bounds), `PeriodicAll` wraps x/y/z.
  pub fn neighbour_atoms_periodic_fixed_positions(&self, cell_id: usize) -> Vec<FixedPositionParticle> {
    let coords = get_coordinates_from_simulation_box_id(cell_id, &self.config.box_count_dim);
    let cell_count = self.config.box_count_dim;
    let world_size = self.config.world_size;
    let (wrap_x, wrap_y, wrap_z) = wrap_flags(self.edge_condition);

    let mut result = Vec::new();

    for x_offset in -1..=1isize {
      let Some((nx, x_shift)) = wrapped_axis(coords.x, x_offset, cell_count.x, wrap_x) else {
        continue;
      };
      for y_offset in -1..=1isize {
        let Some((ny, y_shift)) = wrapped_axis(coords.y, y_offset, cell_count.y, wrap_y) else {
          continue;
        };
        for z_offset in -1..=1isize {
          if x_offset == 0 && y_offset == 0 && z_offset == 0 {
            continue;
          }
          let Some((nz, z_shift)) = wrapped_axis(coords.z, z_offset, cell_count.z, wrap_z) else {
            continue;
          };

          let neighbour_id = get_id_simulation_box(&Vector3::new(nx, ny, nz), &cell_count);
          let shift = Vector3::new(x_shift * world_size.x, y_shift * world_size.y, z_shift * world_size.z);

          for id in self.ids_in_cell(neighbour_id) {
            result.push(FixedPositionParticle { id, position: self.particles[id].get_position() + shift });
          }
        }
      }
    }

    result
  }
}
