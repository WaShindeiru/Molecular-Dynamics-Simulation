use std::collections::HashMap;
use nalgebra::Vector3;
use crate::sim_core::world::boxed_world::box_container::sim_box::{get_coordinates_from_simulation_box_id, get_id_simulation_box};
use crate::sim_core::world::boxed_world::cube::Cube;
use crate::particle::Particle;
use crate::data::constants::get_box_size;
use crate::data::InteractionType;
use crate::output::{AtomDTO, BoxContainerDTO};
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;

pub mod sim_box;
mod integration_cache;

pub struct BoxContainer {
  box_type: InteractionType,
  boxes: Vec<Cube<SimulationBox>>,
  thermostat_epsilon: Vec<f64>,
  box_id_caches: Vec<HashMap<usize, usize>>,
  current_index: usize,

  container_size: Vector3<f64>,
  box_length: Vector3<f64>,

  box_count: usize,
  box_count_dim: Vector3<usize>,

  max_iteration_till_reset: usize,

  integration_box_cache: Cube<SimulationBox>,
  integration_box_id_cache: HashMap<usize, usize>,
  integration_half_velocity_cache: HashMap<usize, Vector3<f64>>,
}

impl BoxContainer {
  pub fn new(atoms: Vec<Particle>, size: Vector3<f64>, box_type: InteractionType,
             max_iteration_till_reset: usize) -> Self {
    let box_length_ = get_box_size(&box_type);
    let box_count_x = (size.x / box_length_).floor();
    let box_count_y = (size.y / box_length_).floor();
    let box_count_z = (size.z / box_length_).floor();

    let box_count_dim = Vector3::new(box_count_x as usize, box_count_y as usize, box_count_z as usize);
    let box_count = box_count_dim.x * box_count_dim.y * box_count_dim.z;

    let box_length_x = size.x / box_count_x;
    let box_length_y = size.y / box_count_y;
    let box_length_z = size.z / box_count_z;
    let box_length = Vector3::new(box_length_x, box_length_y, box_length_z);

    let mut boxes_1: Cube<SimulationBox> = Cube::new(box_count_x as usize, box_count_y as usize, box_count_z as usize);

    for x_i in 0..box_count_x as usize {
      for y_i in 0..box_count_y as usize {
        for z_i in 0..box_count_z as usize {
          let coordinates = Vector3::new(x_i, y_i, z_i);
          let box_id = get_id_simulation_box(&coordinates, &box_count_dim);

          let leftmost_point = Vector3::new(
            x_i as f64 * box_length_x,
            y_i as f64 * box_length_y,
            z_i as f64 * box_length_z,
          );

          let rightmost_point = Vector3::new(
            (x_i + 1) as f64 * box_length_x,
            (y_i + 1) as f64 * box_length_y,
            (z_i + 1) as f64 * box_length_z,
          );

          let sim_box = SimulationBox::new(box_id, leftmost_point, rightmost_point, box_length);

          boxes_1.set(x_i, y_i, z_i, sim_box).unwrap();
        }
      }
    }

    let mut thermostat_epsilon: Vec<f64> = Vec::with_capacity(max_iteration_till_reset);
    thermostat_epsilon.push(0.);

    let mut result = BoxContainer {
      box_type,
      boxes: vec![boxes_1],
      box_id_caches: Vec::new(),
      thermostat_epsilon,
      current_index: 0,

      container_size: size,
      box_length,

      box_count,
      box_count_dim,

      max_iteration_till_reset,

      integration_box_cache: Cube::new(1, 1, 1),
      integration_box_id_cache: HashMap::new(),
      integration_half_velocity_cache: HashMap::new(),
    };

    let box_id_cache = result.assign_particles_to_boxes(atoms);

    let mut box_id_caches: Vec<HashMap<usize, usize>> = Vec::with_capacity(max_iteration_till_reset);
    box_id_caches.push(box_id_cache);

    result.box_id_caches = box_id_caches;

    result
  }

  pub fn assign_box_for_particle(&self, particle: &Particle) -> usize {
    let x: usize;
    let y: usize;
    let z: usize;
    let position = particle.get_position();

    if position.x == self.container_size.x {
      x = self.box_count_dim.x - 1;
    } else {
      x = (position.x / self.box_length.x) as usize;
    }

    if position.y == self.container_size.y {
      y = self.box_count_dim.y - 1;
    } else {
      y = (position.y / self.box_length.y) as usize;
    }

    if position.z == self.container_size.z {
      z = self.box_count_dim.z - 1;
    } else {
      z = (position.z / self.box_length.z) as usize;
    }

    get_id_simulation_box(&Vector3::new(x, y, z), &self.box_count_dim)
  }

  pub fn assign_particles_to_boxes(&mut self, atoms: Vec<Particle>) -> HashMap<usize, usize> {
    let mut box_id_cache: HashMap<usize, usize> = HashMap::with_capacity(atoms.len());

    for particle_i in atoms {

      let box_id = self.assign_box_for_particle(&particle_i);
      box_id_cache.insert(particle_i.get_id() as usize, box_id);
      let coordinates = get_coordinates_from_simulation_box_id(box_id,
                                                               &self.box_count_dim);
      self.boxes.last_mut().unwrap().get_mut(coordinates.x, coordinates.y, coordinates.z)
        .unwrap().add_particle(particle_i);
    }

    box_id_cache
  }

  pub fn set_integration_box_cache(&mut self, cache: HashMap<usize, Particle>) {
    let mut box_id_cache: HashMap<usize, usize> = HashMap::with_capacity(cache.len());
    let mut new_boxes = self.boxes.last().unwrap().clone();
    
    for sim_box in new_boxes.iter_mut() {
      sim_box.clear_box();
    }
    
    for (i_id, particle_i) in cache {
      assert_eq!(i_id, particle_i.get_id() as usize);
      let box_id = self.assign_box_for_particle(&particle_i);
      box_id_cache.insert(particle_i.get_id() as usize, box_id);
      let coordinates = get_coordinates_from_simulation_box_id(box_id,
                                                               &self.box_count_dim);
      new_boxes.get_mut(coordinates.x, coordinates.y, coordinates.z).unwrap()
        .add_particle(particle_i);
    }
    
    self.integration_box_cache = new_boxes;
    self.integration_box_id_cache = box_id_cache;
  }

  pub fn num_of_atoms(&self) -> usize {
    self.box_id_caches.len()
  }

  pub fn box_type(&self) -> &InteractionType {
    &self.box_type
  }

  pub fn container_size(&self) -> &Vector3<f64> {
    &self.container_size
  }

  pub fn box_count(&self) -> usize {
    self.box_count
  }

  pub fn box_count_dim(&self) -> &Vector3<usize> {
    &self.box_count_dim
  }

  pub fn box_length(&self) -> &Vector3<f64> { &self.box_length }

  pub fn current_box_of_id(&self, box_id: usize) -> &SimulationBox {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    self.boxes.last().unwrap().get(coordinates.x, coordinates.y, coordinates.z).unwrap()
  }

  pub fn current_atoms_of_box(&self, box_id: usize) -> &HashMap<usize, Particle> {
    self.current_box_of_id(box_id).particles()
  }

  pub fn current_box_of_atom(&self, atom_id: usize) -> &SimulationBox {
    let box_id = self.box_id_caches.last().unwrap().get(&atom_id).unwrap();
    let coords = get_coordinates_from_simulation_box_id(*box_id, &self.box_count_dim);
    self.boxes.last().unwrap().get(coords.x, coords.y, coords.z).unwrap()
  }

  pub fn box_of_atom_given_index(&self, atom_id: usize, index: usize) -> &SimulationBox {
    let box_id = self.box_id_caches.get(index).unwrap().get(&atom_id).unwrap();
    let coords = get_coordinates_from_simulation_box_id(*box_id, &self.box_count_dim);
    self.boxes.get(index).unwrap().get(coords.x, coords.y, coords.z).unwrap()
  }

  pub fn current_boxes(&self) -> &Cube<SimulationBox> {
    self.boxes.last().unwrap()
  }

  pub fn current_thermostat_epsilon(&self) -> f64 { *self.thermostat_epsilon.last().unwrap() }

  pub fn thermostat_epsilon_of_iteration(&self, iteration: usize) -> f64 {
    *self.thermostat_epsilon.get(iteration).unwrap()
  }
  
  pub fn add_thermostat_epsilon(&mut self, thermostat_epsilon: f64) {
    self.thermostat_epsilon.push(thermostat_epsilon);
  }

  pub fn reset_container(&mut self) {
    let new_index = 0;

    let mut new_thermostat_epsilon: Vec<f64> = Vec::with_capacity(self.max_iteration_till_reset + 1);
    if let Some(last_epsilon) = self.thermostat_epsilon.pop() {
      new_thermostat_epsilon.push(last_epsilon);
    } else {
      panic!("Thermostat epsilon is empty!");
    }

    let mut new_boxes: Vec<Cube<SimulationBox>> = Vec::with_capacity(self.max_iteration_till_reset + 1);
    if let Some(last_boxes) = self.boxes.pop() {
      new_boxes.push(last_boxes);
    } else {
      panic!("Boxes are empty!");
    }

    let mut new_box_id_caches: Vec<HashMap<usize, usize>> = Vec::with_capacity(self.max_iteration_till_reset + 1);
    if let Some(last_box_id_cache) = self.box_id_caches.pop() {
      new_box_id_caches.push(last_box_id_cache);
    } else {
      panic!("Box id caches are empty!");
    }

    self.current_index = new_index;
    self.boxes = new_boxes;
    self.thermostat_epsilon = new_thermostat_epsilon;
    self.box_id_caches = new_box_id_caches;
  }

  pub fn to_transfer_struct(&self, lower_index: usize) -> BoxContainerDTO {
    let mut atoms: Vec<Vec<AtomDTO>> = Vec::new();

    for iteration in lower_index..self.boxes.len() {
      let cube_ = self.boxes.get(iteration).unwrap();

      let mut particles_temp: Vec<AtomDTO> = Vec::new();

      for sim_box in cube_.iter() {
        for (_, particle) in sim_box.particles() {
          particles_temp.push(particle.to_transfer_struct());
        }
      }

      atoms.push(particles_temp);
    }

    BoxContainerDTO {
      box_type: self.box_type,
      atoms,
      thermostat_epsilon: self.thermostat_epsilon.clone(),

      box_length: self.box_length,
      box_count: self.box_count,
      box_count_dim: self.box_count_dim,
    }
  }
}
