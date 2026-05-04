use std::io;
use std::sync::mpsc::{self, TryRecvError};
use std::thread::JoinHandle;

use crate::output::world::boxed::BoxedWorldDTO;
use crate::sim_core::world::saver::PartialWorldSaver;

use super::saver_worker::SaverWorker;

pub(super) struct SaverHandle {
  dto_sender: Option<mpsc::SyncSender<BoxedWorldDTO>>,
  result_receiver: mpsc::Receiver<io::Result<()>>,
  join: Option<JoinHandle<()>>,
}

impl SaverHandle {
  pub(super) fn spawn(saver: PartialWorldSaver, save_enabled: bool) -> Self {
    let (dto_tx, dto_rx) = mpsc::sync_channel::<BoxedWorldDTO>(1);
    let (result_tx, result_rx) = mpsc::channel::<io::Result<()>>();
    let join = if save_enabled {
      Some(std::thread::spawn(move || {
        SaverWorker::new(saver, dto_rx, result_tx).run();
      }))
    } else {
      None
    };

    Self {
      dto_sender: Some(dto_tx),
      result_receiver: result_rx,
      join,
    }
  }

  pub(super) fn shutdown_inner(&mut self) {
    self.dto_sender.take();
    if let Some(j) = self.join.take() {
      let _ = j.join();
    }
  }

  pub(super) fn drain_finished_results(&self) -> io::Result<()> {
    loop {
      match self.result_receiver.try_recv() {
        Ok(Ok(())) => {}
        Ok(Err(e)) => return Err(e),
        Err(TryRecvError::Empty | TryRecvError::Disconnected) => return Ok(()),
      }
    }
  }

  pub(super) fn send_dto(&self, dto: BoxedWorldDTO) -> io::Result<()> {
    let sender = self
      .dto_sender
      .as_ref()
      .ok_or_else(|| io::Error::new(io::ErrorKind::NotConnected, "saver disconnected"))?;

    sender
      .send(dto)
      .map_err(|_| io::Error::new(io::ErrorKind::BrokenPipe, "saver worker terminated"))
  }

  /// Used for final snapshots: waits until the saver has finished persist for this dto.
  pub(super) fn send_dto_wait_result(&self, dto: BoxedWorldDTO) -> io::Result<()> {
    self.send_dto(dto)?;
    match self.result_receiver.recv() {
      Ok(Ok(())) => Ok(()),
      Ok(Err(e)) => Err(e),
      Err(_) => Err(io::Error::new(
        io::ErrorKind::BrokenPipe,
        "saver worker dropped result channel before replying",
      )),
    }
  }
}
