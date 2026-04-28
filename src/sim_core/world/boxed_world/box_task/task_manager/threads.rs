use std::sync::{mpsc, Arc, Mutex};
use std::sync::mpsc::{Receiver, Sender};
use std::thread;
use std::thread::JoinHandle;
use log::info;
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};
use crate::sim_core::world::boxed_world::box_task::handle_task::handle_velocity_batch_task;
use crate::sim_core::world::boxed_world::box_task::task_manager::TaskManager;

pub fn worker_task_handle(job_rx_clone: Arc<Mutex<Receiver<BoxTask>>>, result_tx_clone: Sender<BoxResult>) {
  loop {
    let job = {
      let receiver = job_rx_clone.lock().unwrap();
      receiver.recv()
    };

    match job {
      Ok(task) => {
        match task {
          BoxTask::VelocityBatchTask {
            task_id, box_ids, history, time_step,
            previous_thermostat_epsilon, current_iteration, container_size, edge_condition
          } => {
            let velocity_result = handle_velocity_batch_task(
              task_id, &box_ids, &history, time_step,
              previous_thermostat_epsilon, current_iteration, container_size, edge_condition,
            );
            result_tx_clone.send(BoxResult::VelocityResult(velocity_result)).unwrap();
          },
          
          BoxTask::ForceBatchTask { .. } => {
            unimplemented!("ForceBatchTask not yet implemented");
          }
        }
      }
      Err(_) => {
        break;
      }
    }
  }
}

pub fn create_threads(debug: bool) -> (Sender<BoxTask>, Receiver<BoxResult>, Vec<JoinHandle<()>>, usize) {
  let num_workers: usize;
  if debug {
    num_workers = 1;
  } else  {
    num_workers = thread::available_parallelism()
      .map(|n| n.get() - 1)
      .unwrap_or(1);
  }
  info!("Using {num_workers} threads.");

  let (tx_task, rx_task) = mpsc::channel::<BoxTask>();
  let job_rx = Arc::new(Mutex::new(rx_task));

  let (result_tx, result_rx) = mpsc::channel::<BoxResult>();

  let mut threads: Vec<JoinHandle<()>> = Vec::with_capacity(num_workers);

  for _ in 0..num_workers {
    let job_rx_clone = Arc::clone(&job_rx);
    let result_tx_clone = result_tx.clone();

    let worker = thread::spawn(move || worker_task_handle(job_rx_clone, result_tx_clone));

    threads.push(worker);
  }
  
  (tx_task, result_rx, threads, num_workers)
}

