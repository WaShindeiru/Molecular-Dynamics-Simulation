use std::io;

use crate::persistence::dto::world::boxed::BoxedWorldDTOWithoutHistory;
use crate::sim_core::world::boxed_world::history_manager::HistoryManager;
use crate::sim_core::world::saver::PartialWorldSaver;

pub(crate) type SaverMessage = (BoxedWorldDTOWithoutHistory, HistoryManager);

pub(crate) struct SaverWorker {
  saver: PartialWorldSaver,
  receiver: std::sync::mpsc::Receiver<SaverMessage>,
  result_sender: std::sync::mpsc::Sender<io::Result<()>>,
}

impl SaverWorker {
  pub(crate) fn new(
    saver: PartialWorldSaver,
    receiver: std::sync::mpsc::Receiver<SaverMessage>,
    result_sender: std::sync::mpsc::Sender<io::Result<()>>,
  ) -> Self {
    SaverWorker {
      saver,
      receiver,
      result_sender,
    }
  }

  pub(crate) fn run(mut self) {
    while let Ok((partial, history_manager)) = self.receiver.recv() {
      let lower_index = if partial.number_of_resets > 0 { 1 } else { 0 };
      let dto = history_manager.to_dto(partial, lower_index);

      let atoms = dto
        .history
        .box_container
        .iter()
        .flat_map(|bc| bc.atoms.iter().cloned());
      self.saver.handle_velocity_particles(atoms);

      let result = self.saver.persist_boxed_world(&dto);
      let _ = self.result_sender.send(result);
    }

    let _ = self.saver.persist_velocity_particles();
  }
}
