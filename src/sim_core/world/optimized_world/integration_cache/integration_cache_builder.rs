use std::collections::HashMap;
use std::sync::Arc;
use std::sync::mpsc::{self, Receiver, Sender};
use std::thread;
use std::thread::JoinHandle;

use nalgebra::Vector3;

use crate::sim_core::world::boundary_constraint::{Compliance, ParticleCompliance};
use crate::sim_core::world::boxed_world::box_task::VelocityTaskParticleData;
use crate::sim_core::world::cell::LinkedCellContainer;
use crate::sim_core::world::optimized_world::integration_cache::IntegrationCache;

/// Never read before being overwritten: `add_velocity_results` always writes the real value
/// for id `i` in the same loop iteration where the backing container's `present[i]` becomes
/// `true`, and `build()` only exposes these vectors once both containers report `is_complete`.
fn placeholder_compliance() -> ParticleCompliance {
  ParticleCompliance {
    compliant: false,
    x: Compliance::Compliant,
    y: Compliance::Compliant,
    z: Compliance::Compliant,
  }
}

type PreparedState = (LinkedCellContainer, LinkedCellContainer, Vec<Vector3<f64>>, Vec<ParticleCompliance>);

fn prepare_state(template: &LinkedCellContainer, num_particles: usize) -> PreparedState {
  let local = template.clone();
  let read = template.clone();
  let half_velocity = vec![Vector3::zeros(); num_particles];
  let particle_compliance = vec![placeholder_compliance(); num_particles];
  (local, read, half_velocity, particle_compliance)
}

/// Runs on a single long-lived thread (spawned once in `IntegrationCacheBuilder::new`) instead
/// of being spawned fresh per iteration, so refilling next iteration's scratch state never pays
/// OS thread-creation cost: `build()` just sends a trigger, this loop does the clone, and
/// `ensure_state_ready()` receives the result whenever it's actually needed. Exits cleanly once
/// the builder is dropped and `request_rx` disconnects.
fn prefetch_loop(
  template: Arc<LinkedCellContainer>,
  num_particles: usize,
  request_rx: Receiver<()>,
  response_tx: Sender<PreparedState>,
) {
  while request_rx.recv().is_ok() {
    let state = prepare_state(&template, num_particles);
    if response_tx.send(state).is_err() {
      break;
    }
  }
}

/// Builds the per-iteration `IntegrationCache` (a `local_container` and `read_container`
/// carrying this iteration's post-half-step particle state) while amortizing the cost of
/// allocating those containers' backing storage across the whole simulation run.
///
/// The template cloned by `prefetch_loop` is computed once, via `LinkedCellContainer::reset_clone`,
/// from the very first container the simulation ever sees. It only carries per-particle metadata
/// that is invariant for the lifetime of a run (id, type, mass, ...) — everything
/// iteration-specific is stripped to defaults by `reset_clone`. Because of that invariance, a
/// fresh pair of scratch containers (and fresh `half_velocity`/`particle_compliance` vectors) for
/// iteration N+1 can be produced independent of what iteration N's actual results turn out to
/// be — so all four are prepared together on a long-lived background thread, triggered as soon
/// as the current iteration's state is handed off in `build()`, overlapping with `force_step`'s
/// worker-thread computation instead of sitting on the critical path. The thread is spawned once
/// (not per iteration) so refilling never pays OS thread-creation cost.
pub struct IntegrationCacheBuilder {
  local_container: Option<LinkedCellContainer>,
  read_container: Option<LinkedCellContainer>,
  half_velocity: Vec<Vector3<f64>>,
  particle_compliance: Vec<ParticleCompliance>,
  prefetch_request_tx: Sender<()>,
  prefetch_response_rx: Receiver<PreparedState>,
  _prefetch_thread: JoinHandle<()>,
}

impl IntegrationCacheBuilder {
  pub fn new(old: &LinkedCellContainer) -> Self {
    let num_particles = old.particles().len();
    let template = Arc::new(old.reset_clone());
    let (local_container, read_container, half_velocity, particle_compliance) =
      prepare_state(&template, num_particles);

    let (prefetch_request_tx, request_rx) = mpsc::channel();
    let (response_tx, prefetch_response_rx) = mpsc::channel();
    let _prefetch_thread =
      thread::spawn(move || prefetch_loop(template, num_particles, request_rx, response_tx));

    IntegrationCacheBuilder {
      local_container: Some(local_container),
      read_container: Some(read_container),
      half_velocity,
      particle_compliance,
      prefetch_request_tx,
      prefetch_response_rx,
      _prefetch_thread,
    }
  }

  /// Blocks only if the background prefetch triggered by the previous `build()` call hasn't
  /// finished yet; otherwise this is a cheap `Option` check.
  pub fn ensure_state_ready(&mut self) {
    if self.local_container.is_none() {
      let (local, read, half_velocity, particle_compliance) = self
        .prefetch_response_rx
        .recv()
        .expect("prefetch thread disconnected");
      self.local_container = Some(local);
      self.read_container = Some(read);
      self.half_velocity = half_velocity;
      self.particle_compliance = particle_compliance;
    }
  }

  pub fn add_velocity_results(&mut self, results: HashMap<usize, VelocityTaskParticleData>) {
    let local = self.local_container.as_mut().unwrap();
    let read = self.read_container.as_mut().unwrap();

    for (id, data) in results {
      self.half_velocity[id] = data.half_velocity;
      self.particle_compliance[id] = data.compliance;

      local.change_position(id, data.new_position);
      local.particle_mut(id).set_thermostat_work(data.thermostat_work);

      read.change_position(id, data.new_position);
      read.particle_mut(id).set_thermostat_work(data.thermostat_work);
    }
  }

  /// Replaces converting to a separate dense representation: `is_complete` is a single O(n)
  /// scan over the containers' own presence flags, with no extra allocation or particle copy.
  pub fn build(&mut self) -> Option<IntegrationCache> {
    let ready = self.local_container.as_ref().is_some_and(|c| c.is_complete())
      && self.read_container.as_ref().is_some_and(|c| c.is_complete());
    if !ready {
      return None;
    }

    let local = self.local_container.take().unwrap();
    let read = self.read_container.take().unwrap();

    self.prefetch_request_tx.send(()).expect("prefetch thread disconnected");

    let half_velocity = std::mem::take(&mut self.half_velocity);
    let particle_compliance = std::mem::take(&mut self.particle_compliance);

    Some(IntegrationCache::new(local, Arc::new(read), half_velocity, particle_compliance))
  }
}
