mod thermostat;
pub mod computation;

use std::collections::{HashMap, HashSet};
use std::sync::Arc;
use nalgebra::Vector3;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::computation::FP;
use crate::particle::Particle;
use crate::sim_core::world::boundary_constraint::{check_position_constraint, ParticleCompliance};
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};
use crate::sim_core::world::boxed_world::BoxedWorld;
use crate::sim_core::world::boxed_world::integration::verlet_nose_hoover::thermostat::compute_new_thermostat_epsilon;

impl BoxedWorld {
  pub fn update_verlet_nose_hoover(&mut self, time_step: f64, next_iteration: usize,
                                   desired_temperature: f64, q_effective_mass: f64) {
    let previous_thermostat_epsilon = self.box_container.read().unwrap().current_thermostat_epsilon();

    for sim_box in self.box_container.read().unwrap().current_boxes().iter() {
      if !sim_box.empty() {
        let vel_task = BoxTask::VelocityTask {
          box_container: self.box_container.clone(),
          box_id: sim_box.id(),
          time_step,
          previous_thermostat_epsilon,
          current_iteration: self.iteration,
          container_size: self.size,
        };

        self.tx_task.send(vel_task).unwrap();
      }
    }

    let mut half_velocity_cache_all: HashMap<usize, Vector3<f64>> = HashMap::new();
    let mut particles_all: HashMap<usize, Particle> = HashMap::new();
    let mut compliance_cache: HashMap<usize, ParticleCompliance> = HashMap::new();

    while half_velocity_cache_all.len() < self.num_of_atoms {
      let result = self.rx_result.recv().unwrap();

      match result {
        BoxResult::VelocityResult(task_result) => {
          for (id, result) in task_result.particles {
            half_velocity_cache_all.insert(id, result.half_velocity);
            particles_all.insert(id, result.particle);
            compliance_cache.insert(id, result.compliance);
          }
        },
        _ => {
          panic!("wrong result type");
        }
      }
    }

    assert_eq!(self.num_of_atoms, half_velocity_cache_all.len());
    assert_eq!(self.num_of_atoms, particles_all.len());
    assert_eq!(self.num_of_atoms, compliance_cache.len());

    {
      let mut lock = self.box_container.write().unwrap();
      lock.set_integration_half_velocity_cache(half_velocity_cache_all);
      lock.set_integration_box_cache(particles_all);
    }

    // thermostat epsilon
    let new_thermostat_epsilon;
    {
      let lock = self.box_container.read().unwrap();
      new_thermostat_epsilon =
        compute_new_thermostat_epsilon(previous_thermostat_epsilon, lock.integration_half_velocity_cache(),
                                       lock.atoms_of_integration_box(), time_step, q_effective_mass,
                                       desired_temperature)
    }
    self.box_container.write().unwrap().add_thermostat_epsilon(new_thermostat_epsilon);

    let mut expected_boxes: HashSet<usize> = HashSet::new();

    // add forces computation
    for sim_box in self.box_container.read().unwrap().integration_boxes_cache().iter() {
      if !sim_box.empty() {
        let force_task = BoxTask::ForceTask {
          box_container: Arc::clone(&self.box_container),
          box_id: sim_box.id(),
        };

        self.tx_task.send(force_task).unwrap();

        expected_boxes.insert(sim_box.id());
      }
    }

    let mut box_ids: HashSet<usize> = HashSet::new();

    while box_ids.len() < expected_boxes.len() {
      let result = self.rx_result.recv().unwrap();

      match result {
        BoxResult::ForceResult(result) => {
          let box_id = result.box_id;
          assert!(!box_ids.contains(&box_id));
          box_ids.insert(box_id);

          self.update_integration_cache(&result.acceleration, &result.info_boxed.fp);
        },
        _ => panic!("wrong result type!"),
      }
    }

    assert_eq!(box_ids, expected_boxes);

    // TODO: move setting iteration for a particle into earlier step
    self.box_container.write().unwrap().integration_box_set_velocity(time_step, new_thermostat_epsilon,
                                                                     self.iteration + 1,
                                                                     &compliance_cache);
    // idk
    self.box_container.write().unwrap().apply_integration_cache();

    self.iteration += 1;
    assert_eq!(self.iteration, next_iteration);
  }

  pub fn update_integration_cache(&mut self, acceleration: &HashMap<usize, Vector3<f64>>,
                                  fp: &HashMap<usize, FP>) {
    let mut lock = self.box_container.write().unwrap();

    for (i_id_, i_acc) in acceleration {
      let i_id = *i_id_;
      let i_fp = fp.get(i_id_).unwrap();
      let i_pot_energy = i_fp.potential_energy;
      let i_force = i_fp.force;

      lock.update_integration_box_force(i_id, &i_force, i_acc, i_pot_energy);
    }
  }
}