use std::io;

use clap::{Parser, Subcommand};

use self::combine::{CombineCommand, combine_command};
use self::parse::{ParseCommand, parse_command};
use self::run::{RunCommand, run_simulation_command};

pub mod combine;
pub mod parse;
pub mod run;

/// Run molecular dynamics simulation from JSON config files.
#[derive(Debug, Parser)]
#[command(name = "carbon_nanotube", version, about)]
pub struct Cli {
  #[command(subcommand)]
  pub command: Command,
}

#[derive(Debug, Subcommand)]
pub enum Command {
  Run(RunCommand),
  Parse(ParseCommand),
  Combine(CombineCommand),
}

pub fn run() -> io::Result<()> {
  let cli = Cli::parse();

  match cli.command {
    Command::Run(command) => run_simulation_command(command),
    Command::Parse(command) => parse_command(command),
    Command::Combine(command) => combine_command(command),
  }
}
