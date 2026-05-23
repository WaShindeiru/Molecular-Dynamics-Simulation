use std::collections::HashMap;

use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::{EdgeCondition, ParticleCompliance};
use crate::sim_core::world::boundary_constraint::periodic::apply_velocity_constraint_periodic;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::{HalfVelocityResult, verlet_noose_hoover_half_velocity_position};
use crate::sim_core::world::computation::compute_new_velocity;
use nalgebra::Vector3;

pub struct PartialVelocityStepParticle {
  pub particle: Particle,
  pub compliance: ParticleCompliance,
}

impl From<PartialVelocityStepParticle> for Particle {
  fn from(value: PartialVelocityStepParticle) -> Self {
    value.particle
  }
}

impl From<&PartialVelocityStepParticle> for Particle {
  fn from(value: &PartialVelocityStepParticle) -> Self {
    value.particle.clone()
  }
}

pub struct PartialVelocityStep {
  world_size: Vector3<f64>,
  edge_condition: EdgeCondition,
  particles: HashMap<usize, PartialVelocityStepParticle>,
  thermostat_work: HashMap<usize, f64>,
  thermostat_epsilon: f64,
  max_iter: usize,
  current_iter: usize,
  time_step_coef: f64,
}

impl PartialVelocityStep {
  pub fn new(
    world_size: Vector3<f64>,
    edge_condition: EdgeCondition,
    initial_particles: HashMap<usize, Particle>,
    initial_compliance: HashMap<usize, ParticleCompliance>,
    thermostat_epsilon: f64,
    max_iter: usize,
    time_step: f64,
  ) -> Self {
    let thermostat_work = initial_particles.keys().map(|id| (*id, 0.0)).collect();
    let particles = initial_particles
      .into_iter()
      .map(|(id, particle)| {
        let compliance = *initial_compliance.get(&id).unwrap();
        (id, PartialVelocityStepParticle { particle, compliance })
      })
      .collect();
    Self {
      world_size,
      edge_condition,
      particles,
      thermostat_work,
      thermostat_epsilon,
      max_iter,
      current_iter: 0,
      time_step_coef: time_step / max_iter as f64,
    }
  }

  fn do_one_step(&mut self) -> HashMap<usize, Vector3<f64>> {
    debug_assert!(self.current_iter < self.max_iter);
    let previous_particles = &self.particles;
    let atom_count = previous_particles.len();

    let half_velocity_result = verlet_noose_hoover_half_velocity_position(
      previous_particles.values().map(|f| &f.particle),
      self.time_step_coef,
      self.thermostat_epsilon,
      atom_count,
      1,
      &self.world_size,
      self.edge_condition,
    );

    let half_velocity = half_velocity_result.half_velocity;
    let new_position = half_velocity_result.new_position;
    let thermostat_work = half_velocity_result.thermostat_work;
    let compliance = half_velocity_result.compliance;

    for (i, work) in thermostat_work {
      *self.thermostat_work.get_mut(&i).unwrap() += work;
    }

    let corrected_half_velocity: HashMap<usize, Vector3<f64>> = half_velocity
      .iter()
      .map(|(id, vel)| (*id, apply_velocity_constraint_periodic(compliance.get(id).unwrap(), vel)))
      .collect();

    let correct_compliance: HashMap<usize, ParticleCompliance> = compliance
      .iter()
      .map(|(id, comp)| match comp.compliant {
            true => {
              let prev_comp = previous_particles.get(&id).unwrap().compliance;
              if !prev_comp.compliant {
                (*id, prev_comp)
              } else {
                (*id, *comp)
              }
            },
            false => (*id, *comp),
        })
      .collect();

    let new_velocity: HashMap<usize, Vector3<f64>> = corrected_half_velocity
      .iter()
      .map(|(id, half_vel)| 
        (*id, compute_new_velocity(
          half_vel.clone(), 
          *previous_particles.get(id).unwrap().particle.get_acceleration(), 
          self.thermostat_epsilon, 
          self.time_step_coef)))
      .collect();

    let new_particles: HashMap<usize, PartialVelocityStepParticle> = previous_particles
      .iter()
      .map(|(id, partial_velocity_particle)| {
        let mut new_particle = partial_velocity_particle.particle.clone();
        new_particle.set_velocity(new_velocity.get(id).unwrap().clone());
        new_particle.update_position(new_position.get(id).unwrap().clone());
        (*id, PartialVelocityStepParticle {
          compliance: correct_compliance.get(id).unwrap().clone(),
          particle: new_particle,
        })
      } )
      .collect();

    self.current_iter += 1;
    self.particles = new_particles;

    corrected_half_velocity
  }

  pub fn run(&mut self) -> HalfVelocityResult {
    let mut half_velocity: HashMap<usize, Vector3<f64>> = HashMap::new();

    for _ in 0..self.max_iter {
      half_velocity = self.do_one_step();
    }

    debug_assert!(self.current_iter == self.max_iter);

    let final_particles = &self.particles;

    let mut compliance: HashMap<usize, ParticleCompliance> = HashMap::new();
    let mut new_position: HashMap<usize, Vector3<f64>> = HashMap::new();

    for (i, particle) in final_particles {
      compliance.insert(*i, particle.compliance);
      new_position.insert(*i, *particle.particle.get_position());
    }

    HalfVelocityResult {
      half_velocity,
      new_position,
      thermostat_work: self.thermostat_work.clone(),
      compliance,
    }
  }
}