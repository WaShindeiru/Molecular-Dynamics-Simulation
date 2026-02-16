use std::collections::{HashMap, HashSet};
use nalgebra::Vector3;
use crate::sim_core::world::boxed_world::box_container::sim_box::SimulationBox;
use crate::particle::Particle;
use crate::sim_core::world::boxed_world::box_container::BoxContainer;
use crate::sim_core::world::boxed_world::box_container::sim_box::{get_coordinates_from_simulation_box_id, get_id_simulation_box};
use crate::sim_core::world::boxed_world::cube::Cube;

impl BoxContainer {
  pub fn set_integration_half_velocity_cache(&mut self, cache: HashMap<usize, Vector3<f64>>) {
    self.integration_half_velocity_cache = cache;
  }

  pub fn atoms_of_integration_box(&self) -> impl Iterator<Item = &Particle> {
    self.integration_box_cache.iter().flat_map(|sim_box| sim_box.particles().values())
  }

  pub fn atoms_of_given_integration_box(&self, box_id: usize) -> impl Iterator<Item = &Particle> {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    self.integration_box_cache.get(coordinates.x, coordinates.y, coordinates.z).unwrap().particles().values()
  }

  pub fn integration_half_velocity_cache(&self) -> &HashMap<usize, Vector3<f64>> {
    &self.integration_half_velocity_cache
  }

  pub fn integration_box_id_cache(&self) -> &HashMap<usize, usize> {
    &self.integration_box_id_cache
  }

  pub fn integration_boxes_cache(&self) -> &Cube<SimulationBox> {
    &self.integration_box_cache
  }

  pub fn update_integration_box_force(&mut self, particle_id: usize, force: &Vector3<f64>,
                                      acceleration: &Vector3<f64>, pot_energy: f64) {
    let box_id = self.integration_box_id_cache.get(&particle_id).unwrap();
    let coordinates = get_coordinates_from_simulation_box_id(*box_id, self.box_count_dim());
    let sim_box = self.integration_box_cache.get_vec_mut(&coordinates).unwrap();
    let particle = sim_box.particle_mut(particle_id);
    particle.set_force(particle.get_force() + force);
    particle.set_acceleration(particle.get_acceleration() + acceleration);
    particle.set_potential_energy(particle.get_potential_energy() + pot_energy);
  }

  pub fn neighbours_of_box(&self, box_id: usize) -> Vec<usize> {
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    let mut neighbours: Vec<usize> = Vec::new();
    let box_count_dim = self.box_count_dim();

    for x_ in coordinates.x - 1..coordinates.x+2 {
      for y_ in coordinates.y - 1..coordinates.y+2 {
        for z_ in coordinates.z-1 .. coordinates.z+2 {
          if x_ == coordinates.x && y_ == coordinates.y && z_ == coordinates.z {
            continue;
          }
          if x_ >= 0 && x_ < box_count_dim.x && y_ >= 0 && y_ < box_count_dim.y && z_ > 0 &&
            z_ < box_count_dim.z {
            neighbours.push(get_id_simulation_box(&Vector3::new(x_, y_, z_), box_count_dim));
          }
        }
      }
    }

    neighbours
  }

  pub fn integration_particles_of_neighbour_boxes(&self, box_id: usize) -> impl Iterator<Item = &Particle> {    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    let coordinates = get_coordinates_from_simulation_box_id(box_id, &self.box_count_dim);
    let box_count_dim = self.box_count_dim();

    let x = coordinates.x as isize;
    let y = coordinates.y as isize;
    let z = coordinates.z as isize;

    (-1..=1).flat_map(move |x_offset| {
      (-1..=1).flat_map(move |y_offset| {
        (-1..=1).filter_map(move |z_offset| {
          if x_offset == 0 && y_offset == 0 && z_offset == 0 {
            return None;
          }
          let x_ = x + x_offset;
          let y_ = y + y_offset;
          let z_ = z + z_offset;

          if x_ >= 0 && x_ < box_count_dim.x as isize &&
            y_ >= 0 && y_ < box_count_dim.y as isize &&
            z_ >= 0 && z_ < box_count_dim.z as isize {
            self.integration_box_cache.get(x_ as usize, y_ as usize, z_ as usize)
              .map(|sim_box| sim_box.particles().values())
          } else {
            None
          }
        }).flatten()
      })
    })
  }

  pub fn integration_box_set_velocity(&mut self, time_step: f64, thermostat_epsilon: f64, next_iteration: usize) {
    for sim_box in self.integration_box_cache.iter_mut() {
      for (i_id_, particle_i) in sim_box.particles_mut().iter_mut() {
        let numerator: Vector3<f64> = self.integration_half_velocity_cache.get(i_id_).unwrap() +
          0.5 * particle_i.get_acceleration() * time_step;

        let denominator = 1.0 + 0.5 * time_step * thermostat_epsilon;

        let new_velocity: Vector3<f64> = numerator / denominator;

        particle_i.set_velocity(new_velocity);
        particle_i.set_iteration(next_iteration);
      }
    }
  }

  pub fn apply_integration_cache(&mut self) {
    let empty_cube: Cube<SimulationBox> = Cube::new(1, 1, 1);
    let cache = std::mem::replace(&mut self.integration_box_cache, empty_cube);
    self.boxes.push(cache);

    let empty_id_cache: HashMap<usize, usize> = HashMap::new();
    let id_cache = std::mem::replace(&mut self.integration_box_id_cache, empty_id_cache);
    self.box_id_caches.push(id_cache);

    self.current_index += 1;
    assert_eq!(self.boxes.len() - 1, self.current_index);
    assert_eq!(self.box_id_caches.len() - 1, self.current_index);
    assert_eq!(self.thermostat_epsilon.len() - 1, self.current_index);
  }
}