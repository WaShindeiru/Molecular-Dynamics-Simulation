use std::sync::{mpsc, Arc, Mutex};
use std::sync::mpsc::{Receiver, Sender};
use std::thread;
use std::thread::JoinHandle;
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};
use crate::sim_core::world::boxed_world::box_task::handle_task::handle_half_velocity_position_task;

pub fn create_threads() -> (Sender<BoxTask>, Receiver<BoxResult>, Vec<JoinHandle<()>>) {
  let num_workers = thread::available_parallelism()
    .map(|n| n.get() - 1)
    .unwrap_or(1);

  let (tx_task, rx_task) = mpsc::channel::<BoxTask>();
  let job_rx = Arc::new(Mutex::new(rx_task));

  let (result_tx, result_rx) = mpsc::channel::<BoxResult>();

  let mut threads: Vec<JoinHandle<()>> = Vec::with_capacity(num_workers);

  for _ in 0..num_workers {
    let job_rx_clone = Arc::clone(&job_rx);
    let result_tx_clone = result_tx.clone();

    let worker = thread::spawn(move || {
      loop {
        let job = {
          let receiver = job_rx_clone.lock().unwrap();
          receiver.recv()
        };

        match job {
          Ok(task) => {
            match task {
              BoxTask::VelocityTask {
                box_container, box_id, time_step,
                previous_thermostat_epsilon, current_iteration
              } => {
                let velocity_result = handle_half_velocity_position_task(box_container, 
                                                                   box_id, time_step, 
                                                                   previous_thermostat_epsilon, 
                                                                   current_iteration);
              
                let box_result = BoxResult::VelocityResult {
                  half_velocity_cache: velocity_result.half_velocity_cache,
                  new_position_atoms: velocity_result.new_position_atoms,
                };
                result_tx_clone.send(box_result).unwrap();
              },
              BoxTask::ForceTask => unimplemented!("not yet")
            }
          }
          Err(_) => {
            break;
          }
        }
      }
    });

    threads.push(worker);
  }

  (tx_task, result_rx, threads)
}
