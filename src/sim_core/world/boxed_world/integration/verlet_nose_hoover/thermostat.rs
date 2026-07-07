use crate::data::units::K_B;
use crate::particle::Particle;
use nalgebra::Vector3;
use std::collections::HashMap;

fn compute_half_velocity_kinetic_energy<I>(
  half_velocity_cache: &HashMap<usize, Vector3<f64>>,
  new_position_atoms: I,
) -> (f64, usize)
where
  I: IntoIterator,
  I::Item: AsRef<Particle>,
{
  let mut kinetic_energy = 0.;
  let mut num_of_particles = 0;

  for temp_i in new_position_atoms.into_iter() {
    let particle_i = temp_i.as_ref();
    let i_id = particle_i.get_id();
    let velocity = half_velocity_cache.get(&i_id).unwrap();
    let mass = particle_i.get_mass();

    kinetic_energy += mass * velocity.magnitude().powi(2) / 2.0;
    num_of_particles += 1;
  }

  (kinetic_energy, num_of_particles)
}

pub fn compute_new_thermostat_epsilon<I>(
  thermostat_epsilon: f64,
  half_velocity_cache: &HashMap<usize, Vector3<f64>>,
  new_position_atoms: I,
  time_step: f64,
  q_effective_mass: f64,
  desired_temperature: f64,
) -> f64
where
  I: IntoIterator,
  I::Item: AsRef<Particle>,
{
  let mut new_thermostat_epsilon = thermostat_epsilon;

  let (half_velocity_kinetic_energy, num_of_particles) =
    compute_half_velocity_kinetic_energy(half_velocity_cache, new_position_atoms);

  new_thermostat_epsilon += time_step / q_effective_mass
    * (half_velocity_kinetic_energy
      - 3. / 2. * K_B * desired_temperature * num_of_particles as f64);

  new_thermostat_epsilon
}
