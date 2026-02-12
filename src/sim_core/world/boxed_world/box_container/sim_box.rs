use nalgebra::Vector3;

pub fn get_id_simulation_box(coordinates: &Vector3<usize>, box_count_dim: &Vector3<usize>) -> usize {
  coordinates.x + coordinates.y * box_count_dim.x + coordinates.z * box_count_dim.x * box_count_dim.y
}

pub fn get_coordinates_from_simulation_box_id(box_id: usize, box_count_dim: &Vector3<usize>) -> Vector3<usize> {
  let mut temp = box_id;
  let x = temp % box_count_dim.x;
  temp = temp - x;
  temp = temp / box_count_dim.x;
  let y = temp % box_count_dim.y;
  temp = temp - y;
  let z = temp / box_count_dim.y;
  Vector3::new(x, y, z)
}

#[derive(Clone)]
pub struct SimulationBox {
  id: usize,
  leftmost_point: Vector3<f64>,
  rightmost_point: Vector3<f64>,
  length: Vector3<f64>,
  particle_indexes: Vec<usize>,
}

impl Default for SimulationBox {
  fn default() -> Self {
    SimulationBox::new(
      0,
      Vector3::zeros(),
      Vector3::zeros(),
      Vector3::zeros(),
    )
  }
}

impl SimulationBox {
  pub fn new(id: usize, leftmost_point: Vector3<f64>, rightmost_point: Vector3<f64>,
             length: Vector3<f64>) -> Self {
    SimulationBox {
      id,
      leftmost_point,
      rightmost_point,
      length,
      particle_indexes:  Vec::new(),
    }
  }

  pub fn id(&self) -> usize {
    self.id
  }

  pub fn leftmost_point(&self) -> &Vector3<f64> {
    &self.leftmost_point
  }

  pub fn rightmost_point(&self) -> &Vector3<f64> {
    &self.rightmost_point
  }
  
  pub fn particle_indexes(&self) -> &Vec<usize> {
    &self.particle_indexes
  }
  
  pub fn num_of_particles(&self) -> usize {
    self.particle_indexes.len()
  }
  
  pub fn add_particle_index(&mut self, index: usize) {
    self.particle_indexes.push(index);
  }

  pub fn len(&self) -> usize {
    self.particle_indexes.len()
  }

  pub fn empty(&self) -> bool {
    self.particle_indexes.len() == 0
  }
}

#[cfg(test)]
mod tests {
  use super::{get_coordinates_from_simulation_box_id, get_id_simulation_box};
  use nalgebra::Vector3;

  #[test]
  fn test_get_id_origin() {
    let coords = Vector3::new(0, 0, 0);
    let box_count_dim = Vector3::new(4, 3, 2);

    let expected_id = 0;
    let actual_id = get_id_simulation_box(&coords, &box_count_dim);

    assert_eq!(actual_id, expected_id);

    let actual_coords = get_coordinates_from_simulation_box_id(actual_id, &box_count_dim);

    assert_eq!(actual_coords, coords);
  }

  #[test]
  fn test_get_id_nontrivial() {
    let coords = Vector3::new(2, 1, 3);
    let box_count_dim = Vector3::new(4, 5, 6);

    let expected_id = 66;
    let actual_id = get_id_simulation_box(&coords, &box_count_dim);

    let actual_coordinates = get_coordinates_from_simulation_box_id(actual_id, &box_count_dim);

    assert_eq!(actual_id, expected_id);
    assert_eq!(coords, actual_coordinates);
  }

  #[test]
  fn test_coordinates_roundtrip() {
    let box_count_dim = Vector3::new(3, 4, 5);
    let coordinates = Vector3::new(2, 3, 3);

    let expected_id = 47;
    let actual_id = get_id_simulation_box(&coordinates, &box_count_dim);

    assert_eq!(actual_id, expected_id);

    let coords = get_coordinates_from_simulation_box_id(actual_id, &box_count_dim);

    assert_eq!(coords, coordinates);
  }
}
