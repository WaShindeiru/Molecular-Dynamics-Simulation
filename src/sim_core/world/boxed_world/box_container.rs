use nalgebra::Vector3;
use crate::sim_core::world::boxed_world::box_container::sim_box::get_id_simulation_box;
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

  container_size: Vector3<f64>,
  box_length: Vector3<f64>,

  box_count: usize,
  box_count_dim: Vector3<usize>,

  max_iteration_till_reset: usize,
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

    for (i, particle_i) in atom_container.get_atoms().iter().enumerate() {
      assert_eq!(i, particle_i.get_id() as usize);
      let x: usize;
      let y: usize;
      let z: usize;
      let position = particle_i.get_position();

      if position.x == size.x {
        x = box_count_dim.x - 1;
      } else {
        x = (position.x / box_length_x) as usize;
      }

      if position.y == size.y {
        y = box_count_dim.y - 1;
      } else {
        y = (position.y / box_length_y) as usize;
      }

      if position.z == size.z {
        z = box_count_dim.z - 1;
      } else {
        z = (position.z / box_length_z) as usize;
      }


      boxes_1.get_mut(x, y, z).unwrap().add_particle_index(particle_i.get_id() as usize);
    }

    let mut atoms: Vec<SimpleAtomContainer> = Vec::with_capacity(max_iteration_till_reset);
    atoms.push(atom_container);

    let mut boxes: Vec<Cube<SimulationBox>> = Vec::with_capacity(max_iteration_till_reset);
    boxes.push(boxes_1);

    BoxContainer {
      box_type,
      atoms,
      boxes,

      container_size: size,
      box_length,

      box_count,
      box_count_dim,

      max_iteration_till_reset,
    }
  }

  pub fn box_type(&self) -> &InteractionType {
    &self.box_type
  }

  pub fn container_size(&self) -> &Vector3<f64> {
    &self.container_size
  }

  pub fn get_box_count(&self) -> usize {
    self.box_count
  }

  pub fn get_box_count_dim(&self) -> &Vector3<usize> {
    &self.box_count_dim
  }
}
