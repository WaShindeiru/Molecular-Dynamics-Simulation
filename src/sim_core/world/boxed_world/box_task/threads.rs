use std::sync::{mpsc, Arc, Mutex};
use std::sync::mpsc::{Receiver, Sender};
use std::thread;
use std::thread::JoinHandle;
use crate::sim_core::world::boxed_world::box_task::{BoxResult, BoxTask};
use crate::sim_core::world::boxed_world::box_task::handle_task::{handle_force_task, handle_half_velocity_position_task};

pub fn create_threads(debug: bool) -> (Sender<BoxTask>, Receiver<BoxResult>, Vec<JoinHandle<()>>) {
  let num_workers: usize;
  if debug {
    num_workers = 1;
  } else  {
    num_workers = thread::available_parallelism()
      .map(|n| n.get() - 1)
      .unwrap_or(1); 
  }

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
                previous_thermostat_epsilon, current_iteration, container_size
              } => {
                let velocity_result = handle_half_velocity_position_task(box_container, 
                                                                   box_id, time_step, 
                                                                   previous_thermostat_epsilon, 
                                                                   current_iteration, container_size);
              
                let box_result = BoxResult::VelocityResult(velocity_result);
                result_tx_clone.send(box_result).unwrap();
              },
              BoxTask::ForceTask {
                box_container, box_id,
              } => {
                let force_result = handle_force_task(box_container, box_id);
                
                result_tx_clone.send(BoxResult::ForceResult(force_result)).unwrap();
              }
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
