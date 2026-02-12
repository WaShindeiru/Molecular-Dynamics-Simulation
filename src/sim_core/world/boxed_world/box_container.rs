use std::collections::HashMap;
use nalgebra::Vector3;
use crate::sim_core::world::boxed_world::box_container::sim_box::{get_coordinates_from_simulation_box_id, get_id_simulation_box};
use crate::sim_core::world::boxed_world::cube::Cube;
use crate::particle::Particle;
use crate::data::constants::get_box_size;
use crate::particle::SimpleAtomContainer;
use crate::data::InteractionType;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;

pub mod sim_box;

pub struct BoxContainer {
  box_type: InteractionType,
  atoms: Vec<SimpleAtomContainer>,
  boxes: Vec<Cube<SimulationBox>>,
  thermostat_epsilon: Vec<f64>,
  box_id_caches: Vec<HashMap<usize, usize>>,
  current_index: usize,

  container_size: Vector3<f64>,
  box_length: Vector3<f64>,

  box_count: usize,
  box_count_dim: Vector3<usize>,

  max_iteration_till_reset: usize,

  integration_atoms_cache: Vec<Particle>,
  integration_box_id_cache: HashMap<usize, usize>,
}

impl BoxContainer {
  pub fn new(atoms: Vec<Particle>, size: Vector3<f64>, box_type: InteractionType,
             max_iteration_till_reset: usize) -> Self {
    let box_length_ = get_box_size(&box_type);
    let box_count_x = (size.x / box_length_).floor();
    let box_count_y = (size.y / box_length_).floor();
    let box_count_z = (size.z / box_length_).floor();

    let box_count_dim = Vector3::new(box_count_x as usize, box_count_y as usize, box_count_z as usize);
    let box_count = box_count_dim.x + box_count_dim.y + box_count_dim.z;

    let box_length_x = size.x / box_count_x;
    let box_length_y = size.y / box_count_y;
    let box_length_z = size.z / box_count_z;
    let box_length = Vector3::new(box_length_x, box_length_y, box_length_z);

    let atom_container = SimpleAtomContainer::new_from_atoms(atoms);

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

    let mut atoms: Vec<SimpleAtomContainer> = Vec::with_capacity(max_iteration_till_reset);
    atoms.push(atom_container);

    let mut thermostat_epsilon: Vec<f64> = Vec::with_capacity(max_iteration_till_reset);
    thermostat_epsilon.push(0.);

    let mut result = BoxContainer {
      box_type,
      atoms,
      boxes: Vec::new(),
      box_id_caches: Vec::new(),
      thermostat_epsilon,
      current_index: 0,

      container_size: size,
      box_length,

      box_count,
      box_count_dim,

      max_iteration_till_reset,

      integration_atoms_cache: Vec::new(),
      integration_box_id_cache: HashMap::new(),
    };

    let box_id_cache = result.assign_atoms_to_boxes(&result.atoms.last().unwrap().get_atoms());

    for (particle_id, box_id) in &box_id_cache {
      let coordinates = get_coordinates_from_simulation_box_id(*box_id, &result.box_count_dim);
      boxes_1.get_mut(coordinates.x, coordinates.y, coordinates.z).unwrap().add_particle_index(*particle_id);
    }


    let mut boxes: Vec<Cube<SimulationBox>> = Vec::with_capacity(max_iteration_till_reset);
    boxes.push(boxes_1);

    let mut box_id_caches: Vec<HashMap<usize, usize>> = Vec::with_capacity(max_iteration_till_reset);
    box_id_caches.push(box_id_cache);

    result.boxes = boxes;
    result.box_id_caches = box_id_caches;

    result
  }

  pub fn assign_atoms_to_boxes(&self, atoms: &Vec<Particle>) -> HashMap<usize, usize> {
    let mut box_id_cache: HashMap<usize, usize> = HashMap::with_capacity(atoms.len());

    for (i, particle_i) in atoms.iter().enumerate() {
      assert_eq!(i, particle_i.get_id() as usize);
      let x: usize;
      let y: usize;
      let z: usize;
      let position = particle_i.get_position();

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

      box_id_cache.insert(particle_i.get_id() as usize, get_id_simulation_box(&Vector3::new(x, y, z), &self.box_count_dim));
    }

    box_id_cache
  }

  pub fn num_of_atoms(&self) -> usize {
    self.atoms.len()
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

  pub fn current_atoms(&self) -> &SimpleAtomContainer {
    assert_eq!(self.current_index, self.atoms.len() - 1);
    &self.atoms.last().unwrap()
  }

  pub fn atoms_given_index(&self, index: usize) -> Option<&SimpleAtomContainer> {
    self.atoms.get(index)
  }

  pub fn current_box_of_id(&self, box_id: usize) -> &SimulationBox {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    self.boxes.last().unwrap().get(coordinates.x, coordinates.y, coordinates.z).unwrap()
  }

  pub fn current_atoms_of_box(&self, box_id: usize) -> Vec<Particle> {
    let mut result: Vec<Particle> = Vec::new();
    let atom_indexes = self.current_box_of_id(box_id).particle_indexes();

    for particle_index in atom_indexes {
      result.push(self.atoms.last().unwrap().get_atom(*particle_index).unwrap().clone());
    }

    result
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

  pub fn reset_container(&mut self) {
    let new_index = 0;

    let mut new_atoms: Vec<SimpleAtomContainer> = Vec::with_capacity(self.max_iteration_till_reset + 1);
    if let Some(last_atom_container) = self.atoms.pop() {
      new_atoms.push(last_atom_container);
    } else {
      panic!("Box container is empty!");
    }

    let mut new_thermostat_epsilon: Vec<f64> = Vec::with_capacity(self.max_iteration_till_reset + 1);
    if let Some(last_epsilon) = self.thermostat_epsilon.pop() {
      new_thermostat_epsilon.push(last_epsilon);
    } else {
      panic!("Thermostat epsilon is empty!");
    }

    let mut new_boxes: Vec<Cube<SimulationBox>> = Vec::with_capacity(self.max_iteration_till_reset);
    if let Some(last_boxes) = self.boxes.pop() {
      new_boxes.push(last_boxes);
    } else {
      panic!("Boxes are empty!");
    }

    self.current_index = new_index;
    self.atoms = new_atoms;
    self.boxes = new_boxes;
    self.thermostat_epsilon = new_thermostat_epsilon;
  }
}
