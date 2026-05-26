use crate::data::InteractionType;
use crate::data::constants::get_box_size;
use crate::data::types::AtomType;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::sim_box::get_id_simulation_box;
use nalgebra::Vector3;

pub type SimulationBoxType = InteractionType;

#[derive(Clone, Copy)]
pub struct BoxContainerConfig {
  pub box_type: SimulationBoxType,
  pub box_length: Vector3<f64>,
  pub box_count: usize,
  pub box_count_dim: Vector3<usize>,
  pub world_size: Vector3<f64>,
}

impl BoxContainerConfig {
  pub fn box_coordinates_for_position(&self, position: &Vector3<f64>) -> Vector3<usize> {
    debug_assert!(
      0.0 <= position.x && position.x <= self.world_size.x,
      "position.x={} must be within [0, {}]",
      position.x,
      self.world_size.x
    );
    debug_assert!(
      0.0 <= position.y && position.y <= self.world_size.y,
      "position.y={} must be within [0, {}]",
      position.y,
      self.world_size.y
    );
    debug_assert!(
      0.0 <= position.z && position.z <= self.world_size.z,
      "position.z={} must be within [0, {}]",
      position.z,
      self.world_size.z
    );

    let x = if position.x == self.world_size.x {
      self.box_count_dim.x - 1
    } else {
      (position.x / self.box_length.x) as usize
    };

    let y = if position.y == self.world_size.y {
      self.box_count_dim.y - 1
    } else {
      (position.y / self.box_length.y) as usize
    };

    let z = if position.z == self.world_size.z {
      self.box_count_dim.z - 1
    } else {
      (position.z / self.box_length.z) as usize
    };

    Vector3::new(x, y, z)
  }

  pub fn box_id_for_position(&self, position: &Vector3<f64>) -> usize {
    get_id_simulation_box(
      &self.box_coordinates_for_position(position),
      &self.box_count_dim,
    )
  }
}

pub fn detect_box_type(atoms: &[Particle]) -> SimulationBoxType {
  let mut fe = false;
  let mut c = false;

  for particle in atoms {
    if matches!(particle.get_type(), AtomType::C | AtomType::C_nanotube) {
      c = true;
    } else {
      fe = true;
    }
  }

  if fe && c {
    InteractionType::FeC
  } else if fe {
    InteractionType::FeFe
  } else {
    InteractionType::CC
  }
}

pub fn new_config(atoms: &[Particle], world_size: Vector3<f64>) -> BoxContainerConfig {
  let box_type = detect_box_type(atoms);

  let box_length_ = get_box_size(&box_type);
  let box_count_x = (world_size.x / box_length_).floor();
  let box_count_y = (world_size.y / box_length_).floor();
  let box_count_z = (world_size.z / box_length_).floor();

  let box_count_dim = Vector3::new(
    box_count_x as usize,
    box_count_y as usize,
    box_count_z as usize,
  );
  let box_count = box_count_dim.x * box_count_dim.y * box_count_dim.z;

  let box_length_x = world_size.x / box_count_x;
  let box_length_y = world_size.y / box_count_y;
  let box_length_z = world_size.z / box_count_z;
  let box_length = Vector3::new(box_length_x, box_length_y, box_length_z);

  BoxContainerConfig {
    box_type,
    box_count,
    box_length,
    box_count_dim,
    world_size,
  }
}
