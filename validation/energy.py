import json

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

import argparse
import os

k_b = 8.617333e-5


def load_time_step(path: str) -> float:
  parameters_path = Path(path) / "parameters.json"
  with open(parameters_path, encoding="utf-8") as f:
    parameters = json.load(f)
  return float(parameters["time_step"])


def count_atoms(path: str) -> int:
  particles_initial_path = Path(path) / "particles_initial.json"
  with open(particles_initial_path, encoding="utf-8") as f:
    particles_initial = json.load(f)
  return sum(
    1 for particle in particles_initial["particles"] if particle["particle_type"] == "Atom"
  )


def show_energy_plot(path: str, thermostat: bool, use_time: bool = True, start: int | None = None, end: int | None = None, temp_ylim: float | None = None) -> None:
  energy_data = pd.read_csv(path + '/energy.csv', header=0)
  energy_data = energy_data[energy_data["iteration"] >= 1].reset_index(drop=True)
  if start is not None:
    energy_data = energy_data[energy_data["iteration"] >= start]
  if end is not None:
    energy_data = energy_data[energy_data["iteration"] <= end]
  energy_data = energy_data.reset_index(drop=True)

  iteration = energy_data["iteration"]

  if use_time:
    time_step = load_time_step(path)
    time_elapsed = iteration * time_step
    x_label = "time elapsed [s]"
  else:
    time_elapsed = iteration
    x_label = "iteration"

  kinetic_energy = energy_data["kinetic_energy"]
  potential_energy = energy_data["potential_energy"]
  potential_gravity_energy = energy_data["potential_gravity_energy"]
  total_energy = energy_data["total_energy"]
  if thermostat:
    thermostat_work = energy_data["thermostat_work_total"]
    thermostat_epsilon = energy_data["thermostat_epsilon"]

  if thermostat:
    total_energy_show = thermostat_work + kinetic_energy + potential_energy  + potential_gravity_energy
  else:
    total_energy_show = kinetic_energy + potential_energy + potential_gravity_energy

  plt.plot(time_elapsed, kinetic_energy, label="kinetic energy")
  plt.plot(time_elapsed, potential_energy, label="potential energy")
  plt.plot(time_elapsed, total_energy_show, label="total energy")
  # plt.plot(time_elapsed, total_energy, label="all_energy")
  if thermostat:
    plt.plot(time_elapsed, thermostat_work, label="thermostat work")
  plt.xlabel(x_label)
  plt.ylabel("Energy [eV]")
  # plt.xlim([130000, 140000])
  # plt.ylim(40000, 40100)
  # plt.xlim(500000, 630000)
  # plt.xlim([0, 3e-10])
  plt.title("Energy plot")
  plt.legend()
  plt.savefig(path + '/energy.png')
  plt.show()

  # plt.figure()
  # plt.plot(time_elapsed, kinetic_energy, label="kinetic energy")
  # plt.plot(time_elapsed, potential_energy, label="potential energy")
  # plt.plot(time_elapsed, total_energy_show, label="total energy")
  # plt.xlabel(x_label)
  # plt.ylabel("Energy [eV]")
  # plt.title("Energy plot")
  # plt.legend()
  # plt.savefig(path + '/energy_simple.png')
  # plt.show()

  plt.figure()
  plt.plot(time_elapsed, potential_gravity_energy, label="gravitational pot energy")
  plt.xlabel(x_label)
  plt.ylabel("Energy [eV]")
  plt.title("Gravitational Potential Energy")
  # plt.xlim([0, 3e-10])
  plt.legend()
  plt.savefig(path + '/gravitational_potential_energy.png')
  plt.show()

  if thermostat:
    total_energy_difference = total_energy + thermostat_work - total_energy.iloc[0]
  else:
    total_energy_difference = total_energy - total_energy.iloc[0]

  plt.figure()
  plt.plot(time_elapsed, total_energy_difference, label="Total energy error")

  max_error_label = total_energy_difference.diff().argmax()
  print(max_error_label)
  print(iteration[max_error_label])
  # diff = total_energy_difference[(time_elapsed > 0.7e-11) & (time_elapsed < 0.71e-11)]
  # print(diff)

  plt.xlabel(x_label)
  plt.ylabel("Energy [eV]")
  plt.title("Total energy error")
  # plt.xlim([0, 2.3e-10])
  # plt.ylim([0, 60])
  # plt.xlim([0, 3e-10])
  plt.savefig(path + '/energy_difference.png')
  plt.show()

  # plt.figure()
  # plt.plot(time_elapsed, total_energy_difference, label="Total energy error")

  # plt.xlabel(x_label)
  # plt.ylabel("Energy [eV]")
  # plt.title("Total energy error")
  # plt.xlim([0, 1e-10])
  # plt.ylim([-2, 20])
  # plt.savefig(path + '/energy_difference_small.png')
  # plt.show()

  plt.figure()
  plt.plot(time_elapsed, total_energy_show, label="Total energy")
  plt.xlabel(x_label)
  plt.ylabel("Energy [eV]")
  plt.title("Total energy")
  # plt.xlim([0, 3e-10])
  plt.savefig(path + '/total_energy.png')
  plt.show()

  num_atoms = count_atoms(path)

  if num_atoms > 0:
    mean_kinetic_energy = kinetic_energy / num_atoms
  else:
    mean_kinetic_energy = kinetic_energy  # Fallback if num_atoms not found

  T = 2 / 3 * mean_kinetic_energy / k_b

  plt.figure()
  plt.plot(time_elapsed, T, label="Temperature")
  plt.xlabel(x_label)
  plt.ylabel("Temperature [K]")
  plt.title(f"Temperature")
  # plt.xlim([0, 3e-10])
  if temp_ylim is not None:
    plt.ylim(top=temp_ylim)
  plt.legend()
  plt.savefig(path + '/Temperature.png')
  plt.show()

  if thermostat:
    plt.figure()
    plt.plot(time_elapsed, thermostat_epsilon, label="Thermostat epsilon")
    plt.xlabel(x_label)
    plt.ylabel("thermostat epsilon")
    plt.title("Thermostat epsilon")
    # plt.xlim([0, 3e-10])
    plt.ylim([-1, 1])
    plt.savefig(path + "/thermostat_epsilon.png")
    plt.show()

    # print((T[15000:16000] - 2000).abs().max())
    # print((T[15000:16000] - 2000).abs().min())
    #
    # temp = T[15500:15800]
    # print(((temp - 2000).abs() < 15.0).sum())
    #
    # print((T - 1600).abs().min())


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description="Plot energy data from a simulation output directory.")
  parser.add_argument("path", type=str, help="Path to the simulation output directory")
  parser.add_argument("--use-time", action=argparse.BooleanOptionalAction, default=True, help="Use time as x-axis (default: True)")
  parser.add_argument("--start", type=float, default=0, help="Start iteration (default: 0)")
  parser.add_argument("--end", type=float, default=5e50, help="End iteration (default: 5e50)")
  parser.add_argument("--thermostat", action=argparse.BooleanOptionalAction, default=True, help="Whether thermostat data is present (default: True)")
  parser.add_argument("--temp-ylim", type=float, default=None, help="Upper y-axis limit for temperature plot")
  args = parser.parse_args()

  # output_dir = "../../output"
  # output_dir = "/media/washindeiru/7E442D59442D1585/md"
  # newest_folder = max([os.path.join(output_dir, d) for d in os.listdir(output_dir)], key=os.path.getmtime)

  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/timestamp_investigation/trash/trash_v1_continued"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/timestamp_investigation/triangle/e-17"

  # newest_folder = "../../output/2026-04-14_12-12-07_exp"
  # compare_different_temps("../../output/2026-04-14_12-12-07_exp")

  show_energy_plot(args.path, args.thermostat, use_time=args.use_time, start=args.start, end=args.end, temp_ylim=args.temp_ylim)
  # compare_different_temps(args.path)
