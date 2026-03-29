use flexi_logger::{Duplicate, FileSpec, Logger, LoggerHandle, WriteMode};
use std::sync::OnceLock;
use chrono::{DateTime, Local};

static LOGGER_HANDLE: OnceLock<LoggerHandle> = OnceLock::new();

pub fn get_save_path() -> String {
  let now: DateTime<Local> = Local::now();
  let time_string = now.format("%Y-%m-%d_%H-%M-%S").to_string();
  "../output/".to_string() + &*time_string
}

pub fn init_logging(directory: String) {
  let label = "simulation";

  LOGGER_HANDLE.get_or_init(|| {
    Logger::try_with_env_or_str("info")
      .expect("invalid RUST_LOG setting")
      .log_to_file(FileSpec::default().directory(directory).basename(label))
      .append()
      .write_mode(WriteMode::BufferAndFlush)
      .duplicate_to_stdout(Duplicate::All)
      .format_for_stdout(flexi_logger::colored_default_format)
      .format_for_files(flexi_logger::detailed_format)
      .start()
      .expect("failed to start flexi_logger")
  });
}

