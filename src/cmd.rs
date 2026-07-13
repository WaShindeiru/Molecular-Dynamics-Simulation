use std::io;

use clap::{Parser, Subcommand};

use self::combine::{CombineCommand, combine_command};
use self::convert::{ConvertCommand, convert_command};
use self::inspect::{InspectCommand, inspect_command};
use self::move_particles::{MoveCommand, move_particles_command};
use self::move_particles_minus::{MoveMinusCommand, move_particles_minus_command};
use self::parse::{ParseCommand, parse_command};
use self::reindex::{ReindexCommand, reindex_command};
use self::run::{RunCommand, run_simulation_command};
use self::vel_man_reindex::{VelManReindexCommand, vel_man_reindex_command};

pub mod combine;
pub mod convert;
pub mod inspect;
pub mod move_particles;
pub mod move_particles_minus;
pub mod parse;
pub mod reindex;
pub mod run;
pub mod vel_man_reindex;

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
  Convert(ConvertCommand),
  Move(MoveCommand),
  MoveMinus(MoveMinusCommand),
  Inspect(InspectCommand),
  Reindex(ReindexCommand),
  VelManReindex(VelManReindexCommand),
}

pub fn run() -> io::Result<()> {
  let cli = Cli::parse();

  match cli.command {
    Command::Run(command) => run_simulation_command(command),
    Command::Parse(command) => parse_command(command),
    Command::Combine(command) => combine_command(command),
    Command::Convert(command) => convert_command(command),
    Command::Move(command) => move_particles_command(command),
    Command::MoveMinus(command) => move_particles_minus_command(command),
    Command::Inspect(command) => inspect_command(command),
    Command::Reindex(command) => reindex_command(command),
    Command::VelManReindex(command) => vel_man_reindex_command(command),
  }
}
