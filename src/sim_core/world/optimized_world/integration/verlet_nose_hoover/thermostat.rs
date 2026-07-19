use nalgebra::Vector3;

use crate::data::units::K_B;
use crate::particle::Particle;

/// Same as `crate::sim_core::world::linked_cell_world::integration::verlet_nose_hoover::
/// thermostat::compute_new_thermostat_epsilon`, but takes a plain `&Particle` iterator instead
/// of `I::Item: AsRef<Particle>` — the dense `LinkedCellContainer` yields
/// `&Particle` directly and `Particle` has no `AsRef<Particle>` impl to bridge that generic bound.
pub fn compute_new_thermostat_epsilon<'a>(
  thermostat_epsilon: f64,
  half_velocity_cache: &[Vector3<f64>],
  new_position_atoms: impl Iterator<Item = &'a Particle>,
  time_step: f64,
  q_effective_mass: f64,
  desired_temperature: f64,
) -> f64 {
  let mut num_of_particles = 0;
  let kinetic_energy: f64 = new_position_atoms
    .map(|p| {
      let vel = half_velocity_cache[p.get_id()];
      num_of_particles += 1;
      p.get_mass() * vel.magnitude().powi(2) / 2.0
    })
    .sum();

  thermostat_epsilon
    + time_step / q_effective_mass
      * (kinetic_energy - 3. / 2. * K_B * desired_temperature * num_of_particles as f64)
}
