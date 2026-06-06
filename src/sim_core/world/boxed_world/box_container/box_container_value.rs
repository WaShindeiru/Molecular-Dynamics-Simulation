use std::collections::HashSet;

use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::box_container_config::BoxContainerConfig;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimBoxEdge::{
  LeftEdge, Normal, RightEdge,
};
use crate::sim_core::world::boxed_world::box_container::sim_box::{
  SimBoxPlacement, SimulationBox, get_coordinates_from_simulation_box_id, get_id_simulation_box,
};
use crate::utils::cube::Cube;
use nalgebra::Vector3;
use std::collections::HashMap;
use std::sync::Arc;

impl BoxContainer<SimulationBox> {
  pub fn new_local(config: BoxContainerConfig) -> Self {
    let mut boxes: Cube<SimulationBox> = Cube::new(
      config.box_count_dim.x,
      config.box_count_dim.y,
      config.box_count_dim.z,
    );

    for x_i in 0..config.box_count_dim.x {
      for y_i in 0..config.box_count_dim.y {
        for z_i in 0..config.box_count_dim.z {
          let coordinates = Vector3::new(x_i, y_i, z_i);
          let box_id = get_id_simulation_box(&coordinates, &config.box_count_dim);

          let leftmost_point = Vector3::new(
            x_i as f64 * config.box_length.x,
            y_i as f64 * config.box_length.y,
            z_i as f64 * config.box_length.z,
          );

          let rightmost_point = Vector3::new(
            (x_i + 1) as f64 * config.box_length.x,
            (y_i + 1) as f64 * config.box_length.y,
            (z_i + 1) as f64 * config.box_length.z,
          );

          let x_edge = if x_i == 0 {
            LeftEdge
          } else if x_i == config.box_count_dim.x - 1 {
            RightEdge
          } else {
            Normal
          };

          let y_edge = if y_i == 0 {
            LeftEdge
          } else if y_i == config.box_count_dim.y - 1 {
            RightEdge
          } else {
            Normal
          };

          let z_edge = if z_i == 0 {
            LeftEdge
          } else if z_i == config.box_count_dim.z - 1 {
            RightEdge
          } else {
            Normal
          };

          let sim_box = SimulationBox::new(
            box_id,
            leftmost_point,
            rightmost_point,
            config.box_length,
            SimBoxPlacement {
              x: x_edge,
              y: y_edge,
              z: z_edge,
            },
          );

          boxes.set(x_i, y_i, z_i, sim_box).unwrap();
        }
      }
    }

    BoxContainer {
      config,
      simulation_boxes: boxes,
      box_id_cache: HashMap::new(),
    }
  }

  pub fn new_local_with_particles_with_box_ids<I: IntoArcParticle>(
    config: BoxContainerConfig,
    particles: &[I],
  ) -> (Self, HashSet<usize>) {
    let mut box_ids: HashSet<usize> = HashSet::new();
    
    let mut container = Self::new_local(config);
    for particle in particles {
      let box_id = container.add_particle(particle.into_arc_particle());
      box_ids.insert(box_id);
    }
    (container, box_ids)
  }

  pub fn new_local_with_particles<I: IntoArcParticle>(
    config: BoxContainerConfig,
    particles: &[I],
  ) -> Self {
    Self::new_local_with_particles_with_box_ids(config, particles).0
  }

  pub fn get_box_mut(&mut self, box_id: usize) -> &mut SimulationBox {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
    self
      .simulation_boxes
      .get_mut(coordinates.x, coordinates.y, coordinates.z)
      .unwrap()
  }

  pub fn add_particle(&mut self, particle: Arc<Particle>) -> usize {
    let particle_id = particle.get_id();
    let box_coordinates = self
      .config()
      .box_coordinates_for_position(particle.get_position());
    self
      .simulation_boxes
      .get_mut(box_coordinates.x, box_coordinates.y, box_coordinates.z)
      .unwrap()
      .add_particle(particle);

    let box_id = get_id_simulation_box(&box_coordinates, &self.config().box_count_dim);
    self.box_id_cache.insert(particle_id, box_id);

    box_id
  }

  pub fn add_particle_with_box_id(&mut self, particle: Arc<Particle>, box_id: usize) {
    let particle_id = particle.get_id();

    #[cfg(debug_assertions)]
    {
      let box_coordinates = self
        .config()
        .box_coordinates_for_position(particle.get_position());
      let computed_box_id =
        get_id_simulation_box(&box_coordinates, &self.config().box_count_dim);
      assert_eq!(
        computed_box_id, box_id,
        "particle {particle_id} box mapping mismatch"
      );
    }

    self.get_box_mut(box_id).add_particle(particle);
    self.box_id_cache.insert(particle_id, box_id);
  }

  pub fn get_box(&self, box_id: usize) -> &SimulationBox {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.config.box_count_dim);
    self
      .simulation_boxes
      .get(coordinates.x, coordinates.y, coordinates.z)
      .unwrap()
  }

  pub fn simulation_boxes(&self) -> &Cube<SimulationBox> {
    &self.simulation_boxes
  }

  pub fn simulation_boxes_mut(&mut self) -> &mut Cube<SimulationBox> {
    &mut self.simulation_boxes
  }

  pub fn into_shared(self) -> BoxContainer<Arc<SimulationBox>> {
    let (x, y, z) = self.simulation_boxes.dimensions();
    let mut shared_boxes: Cube<Arc<SimulationBox>> = Cube::new(x, y, z);

    for ((xi, yi, zi), sim_box) in self.simulation_boxes.iter_with_coords() {
      shared_boxes
        .set(xi, yi, zi, Arc::new(sim_box.clone()))
        .unwrap();
    }

    BoxContainer {
      config: self.config,
      simulation_boxes: shared_boxes,
      box_id_cache: self.box_id_cache,
    }
  }
}

pub trait IntoArcParticle {
  fn into_arc_particle(&self) -> Arc<Particle>;
}

impl IntoArcParticle for Particle {
  fn into_arc_particle(&self) -> Arc<Particle> {
    Arc::new(self.clone())
  }
}

impl IntoArcParticle for Arc<Particle> {
  fn into_arc_particle(&self) -> Arc<Particle> {
    Arc::clone(self)
  }
}
