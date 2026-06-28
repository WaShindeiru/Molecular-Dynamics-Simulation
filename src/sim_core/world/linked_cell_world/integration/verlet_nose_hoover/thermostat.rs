use nalgebra::Vector3;

use crate::data::units::K_B;
use crate::particle::Particle;

pub fn compute_new_thermostat_epsilon<I>(
  thermostat_epsilon: f64,
  half_velocity_cache: &[Vector3<f64>],
  new_position_atoms: I,
  time_step: f64,
  q_effective_mass: f64,
  desired_temperature: f64,
) -> f64
where
  I: IntoIterator,
  I::Item: AsRef<Particle>,
{
  let num_of_particles = half_velocity_cache.len();
  let kinetic_energy: f64 = new_position_atoms
    .into_iter()
    .map(|p| {
      let p = p.as_ref();
      let vel = half_velocity_cache[p.get_id()];
      p.get_mass() * vel.magnitude().powi(2) / 2.0
    })
    .sum();

  thermostat_epsilon
    + time_step / q_effective_mass
      * (kinetic_energy - 3. / 2. * K_B * desired_temperature * num_of_particles as f64)
}
