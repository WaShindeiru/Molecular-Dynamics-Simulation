use crate::particle::Particle;
use crate::persistence::dto::world::boxed::box_container::BoxContainerDTO;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimBoxEdge::{
  LeftEdge, Normal, RightEdge,
};
use crate::sim_core::world::boxed_world::box_container::sim_box::{
  SimBoxPlacement, SimulationBox, get_coordinates_from_simulation_box_id, get_id_simulation_box,
};
use crate::sim_core::world::boxed_world::box_container::{BoxContainer, box_container_config};
use crate::utils::cube::Cube;
use nalgebra::Vector3;
use std::collections::{HashMap, HashSet};
use std::sync::Arc;

impl BoxContainer<Arc<SimulationBox>> {
  pub fn new(atoms: Vec<Particle>, world_size: Vector3<f64>) -> Self {
    let config = box_container_config::new_config(&atoms, world_size);
    let mut container = Self::from_config(config);
    container.box_id_cache = container.assign_particles_to_boxes(atoms);
    container
  }

  fn from_config(box_container_config: BoxContainerConfig) -> Self {
    let mut boxes: Cube<Arc<SimulationBox>> = Cube::new(
      box_container_config.box_count_dim.x,
      box_container_config.box_count_dim.y,
      box_container_config.box_count_dim.z,
    );

    for x_i in 0..box_container_config.box_count_dim.x {
      for y_i in 0..box_container_config.box_count_dim.y {
        for z_i in 0..box_container_config.box_count_dim.z {
          let coordinates = Vector3::new(x_i, y_i, z_i);
          let box_id = get_id_simulation_box(&coordinates, &box_container_config.box_count_dim);

          let leftmost_point = Vector3::new(
            x_i as f64 * box_container_config.box_length.x,
            y_i as f64 * box_container_config.box_length.y,
            z_i as f64 * box_container_config.box_length.z,
          );

          let rightmost_point = Vector3::new(
            (x_i + 1) as f64 * box_container_config.box_length.x,
            (y_i + 1) as f64 * box_container_config.box_length.y,
            (z_i + 1) as f64 * box_container_config.box_length.z,
          );

          let x_edge = if x_i == 0 {
            LeftEdge
          } else if x_i == box_container_config.box_count_dim.x - 1 {
            RightEdge
          } else {
            Normal
          };

          let y_edge = if y_i == 0 {
            LeftEdge
          } else if y_i == box_container_config.box_count_dim.y - 1 {
            RightEdge
          } else {
            Normal
          };

          let sim_box = SimulationBox::new(
            box_id,
            leftmost_point,
            rightmost_point,
            box_container_config.box_length,
            SimBoxPlacement {
              x: x_edge,
              y: y_edge,
            },
          );

          boxes.set(x_i, y_i, z_i, Arc::new(sim_box)).unwrap();
        }
      }
    }

    BoxContainer {
      config: box_container_config,
      simulation_boxes: boxes,
      box_id_cache: HashMap::new(),
    }
  }

  pub fn assign_particles_to_boxes(&mut self, atoms: Vec<Particle>) -> HashMap<usize, usize> {
    let mut box_id_cache: HashMap<usize, usize> = HashMap::with_capacity(atoms.len());

    for particle_i in atoms {
      let box_id = self.assign_box_id_for_particle(&particle_i);
      box_id_cache.insert(particle_i.get_id(), box_id);
      let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
      Arc::make_mut(
        self
          .simulation_boxes
          .get_mut(coordinates.x, coordinates.y, coordinates.z)
          .unwrap(),
      )
      .add_particle(Arc::new(particle_i));
    }

    box_id_cache
  }

  pub fn simulation_boxes(&self) -> &Cube<Arc<SimulationBox>> {
    &self.simulation_boxes
  }

  pub fn get_box(&self, box_id: usize) -> Arc<SimulationBox> {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
    self
      .simulation_boxes
      .get(coordinates.x, coordinates.y, coordinates.z)
      .unwrap()
      .clone()
  }

  // TODO: change this so that don't have to iterate over whole cube maybe?
  pub fn view_select_boxes(&self, box_ids: &[usize]) -> BoxContainer<Option<Arc<SimulationBox>>> {
    let id_set: HashSet<usize> = box_ids.iter().copied().collect();

    let (nx, ny, nz) = self.simulation_boxes.dimensions();
    let mut result: Cube<Option<Arc<SimulationBox>>> = Cube::new(nx, ny, nz);

    for ((x, y, z), sim_box) in self.simulation_boxes.iter_with_coords() {
      let box_id = get_id_simulation_box(&Vector3::new(x, y, z), &self.config.box_count_dim);
      let value = if id_set.contains(&box_id) {
        Some(Arc::clone(sim_box))
      } else {
        None
      };
      result.set(x, y, z, value).unwrap();
    }

    BoxContainer {
      config: self.config,
      simulation_boxes: result,
      box_id_cache: self.box_id_cache.clone(),
    }
  }

  // TODO: Consider if it wouldn't be better to do Arc::make_mut instead?
  pub fn from_particles(
    config: BoxContainerConfig,
    particles: &HashMap<usize, Particle>,
    box_mapping: &HashMap<usize, usize>,
  ) -> Self {
    let mut local = BoxContainer::<SimulationBox>::new_local(config);

    for (particle_id, particle) in particles {
      let box_id = *box_mapping
        .get(particle_id)
        .expect("Particle not found in box mapping");
      local
        .get_box_mut(box_id)
        .add_particle(Arc::new(particle.clone()));
    }

    let mut shared = local.into_shared();
    shared.box_id_cache = box_mapping.clone();
    shared
  }

  pub fn to_transfer_struct(&self) -> BoxContainerDTO {
    let atoms = self
      .simulation_boxes
      .iter()
      .flat_map(|sim_box| sim_box.particles().values())
      .map(|particle| particle.to_transfer_struct())
      .collect();

    BoxContainerDTO {
      atoms,
      config: self.config,
    }
  }
}
