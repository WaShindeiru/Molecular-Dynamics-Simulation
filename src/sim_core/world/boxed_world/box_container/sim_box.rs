use nalgebra::Vector3;

pub fn get_id_simulation_box(coordinates: &Vector3<usize>, box_count_dim: &Vector3<usize>) -> usize {
  coordinates.x + coordinates.y * box_count_dim.x + coordinates.z * box_count_dim.x * box_count_dim.y
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
}
