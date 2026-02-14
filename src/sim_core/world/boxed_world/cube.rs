use nalgebra::Vector3;

#[derive(Clone)]
pub struct Cube<T> {
  data: Vec<Vec<Vec<T>>>,
  size: Vector3<usize>,
}

impl<T> Cube<T> {
  /// Creates a new empty Cube with the specified dimensions
  pub fn new(x_size: usize, y_size: usize, z_size: usize) -> Self
  where
    T: Default + Clone,
  {
    let data = vec![vec![vec![T::default(); z_size]; y_size]; x_size];
    Self {
      data,
      size: Vector3::new(x_size, y_size, z_size),
    }
  }

  /// Creates a new Cube with the specified dimensions filled with the given value
  pub fn new_with_value(x_size: usize, y_size: usize, z_size: usize, value: T) -> Self
  where
    T: Clone,
  {
    let data = vec![vec![vec![value.clone(); z_size]; y_size]; x_size];
    Self {
      data,
      size: Vector3::new(x_size, y_size, z_size),
    }
  }

  /// Creates a new Cube from existing 3D Vec data
  pub fn from_vec(data: Vec<Vec<Vec<T>>>) -> Option<Self> {
    if data.is_empty() {
      return None;
    }

    let x_size = data.len();
    let y_size = data[0].len();
    let z_size = if y_size > 0 { data[0][0].len() } else { 0 };

    // Validate that all dimensions are consistent
    for x in &data {
      if x.len() != y_size {
        return None;
      }
      for y in x {
        if y.len() != z_size {
          return None;
        }
      }
    }

    Some(Self {
      data,
      size: Vector3::new(x_size, y_size, z_size),
    })
  }

  /// Gets a reference to the element at the specified coordinates
  pub fn get(&self, x: usize, y: usize, z: usize) -> Option<&T> {
    self.data.get(x)?.get(y)?.get(z)
  }
  
  pub fn get_vec(&self, vec: &Vector3<usize>) -> Option<&T> {
    self.data.get(vec.x)?.get(vec.y)?.get(vec.z)
  }

  pub fn get_vec_mut(&mut self, vec: &Vector3<usize>) -> Option<&mut T> {
    self.data.get_mut(vec.x)?.get_mut(vec.y)?.get_mut(vec.z)
  }

  /// Gets a mutable reference to the element at the specified coordinates
  pub fn get_mut(&mut self, x: usize, y: usize, z: usize) -> Option<&mut T> {
    self.data.get_mut(x)?.get_mut(y)?.get_mut(z)
  }

  /// Sets the value at the specified coordinates
  pub fn set(&mut self, x: usize, y: usize, z: usize, value: T) -> Result<(), &'static str> {
    if x >= self.size.x || y >= self.size.y || z >= self.size.z {
      return Err("Index out of bounds");
    }
    self.data[x][y][z] = value;
    Ok(())
  }

  /// Returns the dimensions of the cube as (x_size, y_size, z_size)
  pub fn dimensions(&self) -> (usize, usize, usize) {
    (self.size.x, self.size.y, self.size.z)
  }

  /// Returns the x dimension size
  pub fn x_size(&self) -> usize {
    self.size.x
  }

  /// Returns the y dimension size
  pub fn y_size(&self) -> usize {
    self.size.y
  }

  /// Returns the z dimension size
  pub fn z_size(&self) -> usize {
    self.size.z
  }

  /// Returns the total number of elements in the cube
  pub fn len(&self) -> usize {
    self.size.x * self.size.y * self.size.z
  }

  /// Returns true if the cube contains no elements
  pub fn is_empty(&self) -> bool {
    self.len() == 0
  }

  /// Returns a reference to the underlying 3D Vec
  pub fn as_vec(&self) -> &Vec<Vec<Vec<T>>> {
    &self.data
  }

  /// Returns a mutable reference to the underlying 3D Vec
  pub fn as_vec_mut(&mut self) -> &mut Vec<Vec<Vec<T>>> {
    &mut self.data
  }

  /// Consumes the Cube and returns the underlying 3D Vec
  pub fn into_vec(self) -> Vec<Vec<Vec<T>>> {
    self.data
  }

  /// Iterates over all elements in the cube (in x, y, z order)
  pub fn iter(&self) -> impl Iterator<Item = &T> {
    self.data.iter().flat_map(|x| x.iter().flat_map(|y| y.iter()))
  }

  /// Iterates over all elements in the cube with their coordinates
  pub fn iter_with_coords(&self) -> impl Iterator<Item = ((usize, usize, usize), &T)> {
    self.data.iter().enumerate().flat_map(|(x, x_vec)| {
      x_vec.iter().enumerate().flat_map(move |(y, y_vec)| {
        y_vec.iter().enumerate().map(move |(z, item)| ((x, y, z), item))
      })
    })
  }

  /// Mutably iterates over all elements in the cube (in x, y, z order)
  pub fn iter_mut(&mut self) -> impl Iterator<Item = &mut T> {
    self.data.iter_mut().flat_map(|x| x.iter_mut().flat_map(|y| y.iter_mut()))
  }

  /// Fills the entire cube with the given value
  pub fn fill(&mut self, value: T)
  where
    T: Clone,
  {
    for x in 0..self.size.x {
      for y in 0..self.size.y {
        for z in 0..self.size.z {
          self.data[x][y][z] = value.clone();
        }
      }
    }
  }

  /// Clears the cube by filling it with default values
  pub fn clear(&mut self)
  where
    T: Default,
  {
    for x in 0..self.size.x {
      for y in 0..self.size.y {
        for z in 0..self.size.z {
          self.data[x][y][z] = T::default();
        }
      }
    }
  }
}