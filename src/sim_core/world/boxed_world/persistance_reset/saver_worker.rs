use std::io;

use crate::output::world::boxed::BoxedWorldDTO;
use crate::sim_core::world::saver::PartialWorldSaver;

pub(crate) struct SaverWorker {
  saver: PartialWorldSaver,
  receiver: std::sync::mpsc::Receiver<BoxedWorldDTO>,
  result_sender: std::sync::mpsc::Sender<io::Result<()>>,
}

impl SaverWorker {
  pub(crate) fn new(
    saver: PartialWorldSaver,
    receiver: std::sync::mpsc::Receiver<BoxedWorldDTO>,
    result_sender: std::sync::mpsc::Sender<io::Result<()>>,
  ) -> Self {
    SaverWorker {
      saver,
      receiver,
      result_sender,
    }
  }

  pub(crate) fn run(mut self) {
    while let Ok(dto) = self.receiver.recv() {
      let result = self.saver.persist_boxed_world(&dto);
      let _ = self.result_sender.send(result);
    }
  }
}
