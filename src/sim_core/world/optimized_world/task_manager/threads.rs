use nalgebra::Vector3;
use std::sync::mpsc::{Receiver, Sender};
use std::sync::{Arc, Mutex, mpsc};
use std::thread;
use std::thread::JoinHandle;

use crate::sim_core::world::computation::FP;
use crate::sim_core::world::optimized_world::handle_task::{
  handle_force_batch_task, handle_velocity_batch_task,
};
use crate::sim_core::world::optimized_world::optimized_task::{OptimizedResult, OptimizedTask};

pub struct Worker {
  pub fp: Vec<FP>,
  pub gradients_cache: Vec<Vector3<f64>>,
}

impl Worker {
  pub fn new(num_atoms: usize) -> Self {
    Worker {
      fp: vec![FP { force: Vector3::zeros(), potential_energy: 0. }; num_atoms],
      gradients_cache: vec![Vector3::zeros(); num_atoms],
    }
  }
}

fn worker_task_handle(
  job_rx: Arc<Mutex<Receiver<OptimizedTask>>>,
  result_tx: Sender<OptimizedResult>,
  num_atoms: usize,
) {
  let mut worker = Worker::new(num_atoms);

  loop {
    let job = {
      let receiver = job_rx.lock().unwrap();
      receiver.recv()
    };

    match job {
      Ok(task) => match task {
        OptimizedTask::VelocityBatchTask {
          task_id,
          cell_ids,
          history,
          time_step,
          previous_thermostat_epsilon,
          current_iteration,
        } => {
          let velocity_result = handle_velocity_batch_task(
            task_id,
            &*cell_ids,
            &history,
            time_step,
            previous_thermostat_epsilon,
            current_iteration,
          );
          result_tx.send(OptimizedResult::VelocityResult(velocity_result)).unwrap();
        }

        OptimizedTask::ForceBatchTask { task_id, cell_ids, integration_cache } => {
          let force_result = handle_force_batch_task(
            task_id,
            &*cell_ids,
            &integration_cache,
            &mut worker.fp,
            &mut worker.gradients_cache,
          );
          result_tx.send(OptimizedResult::ForceResult(force_result)).unwrap();
        }
      },
      Err(_) => break,
    }
  }
}

pub fn create_threads(
  debug: bool,
  num_atoms: usize,
) -> (
  Sender<OptimizedTask>,
  Receiver<OptimizedResult>,
  Vec<JoinHandle<()>>,
  usize,
) {
  let num_workers = if debug {
    1
  } else {
    thread::available_parallelism().map(|n| n.get() - 1).unwrap_or(1)
  };
  log::info!("Using {num_workers} threads.");

  let (tx_task, rx_task) = mpsc::channel::<OptimizedTask>();
  let job_rx = Arc::new(Mutex::new(rx_task));
  let (result_tx, result_rx) = mpsc::channel::<OptimizedResult>();

  let mut threads = Vec::with_capacity(num_workers);
  for _ in 0..num_workers {
    let job_rx_clone = Arc::clone(&job_rx);
    let result_tx_clone = result_tx.clone();
    threads.push(thread::spawn(move || {
      worker_task_handle(job_rx_clone, result_tx_clone, num_atoms)
    }));
  }

  (tx_task, result_rx, threads, num_workers)
}
