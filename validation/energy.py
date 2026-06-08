import json

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

k_b = 8.617333e-5


def load_time_step(path: str) -> float:
  parameters_path = Path(path) / "parameters.json"
  with open(parameters_path, encoding="utf-8") as f:
    parameters = json.load(f)
  return float(parameters["time_step"])


def show_energy_plot(path: str, thermostat: bool, start: int | None = None, end: int | None = None) -> None:
  energy_data = pd.read_csv(path + '/energy.csv', header=0)
  energy_data = energy_data[energy_data["iteration"] >= 1].reset_index(drop=True)
  if start is not None:
    energy_data = energy_data[energy_data["iteration"] >= start]
  if end is not None:
    energy_data = energy_data[energy_data["iteration"] <= end]
  energy_data = energy_data.reset_index(drop=True)

  iteration = energy_data["iteration"]
  time_step = load_time_step(path)
  time_elapsed = iteration * time_step
  kinetic_energy = energy_data["kinetic_energy"]
  potential_energy = energy_data["potential_energy"]
  # phantom_energy = energy_data["phantom_energy"]
  potential_gravity_energy = energy_data["potential_gravity_energy"]
  total_energy = energy_data["total_energy"]
  if thermostat:
    thermostat_work = energy_data["thermostat_work_total"]
    thermostat_epsilon = energy_data["thermostat_epsilon"]

  total_energy_show = thermostat_work + kinetic_energy + potential_energy  + potential_gravity_energy

  plt.plot(time_elapsed, kinetic_energy, label="kinetic energy")
  plt.plot(time_elapsed, potential_energy, label="potential energy")
  # plt.plot(time_elapsed, phantom_energy, label="phantom energy")
  plt.plot(time_elapsed, total_energy_show, label="total energy")
  # plt.plot(time_elapsed, total_energy, label="all_energy")
  if thermostat:
    plt.plot(time_elapsed, thermostat_work, label="thermostat work")
  plt.xlabel("time elapsed [s]")
  plt.ylabel("Energy [eV]")
  # plt.xlim([130000, 140000])
  # plt.ylim(40000, 40100)
  # plt.xlim(500000, 630000)
  plt.title("Energy plot")
  plt.legend()
  plt.savefig(path + '/energy.png')
  plt.show()

  plt.figure()
  plt.plot(time_elapsed, potential_gravity_energy, label="gravitational pot energy")
  plt.xlabel("time elapsed [s]")
  plt.ylabel("Energy [eV]")
  plt.title("Gravitational Potential Energy")
  plt.legend()
  plt.savefig(path + '/gravitational_potential_energy.png')
  plt.show()

  if thermostat:
    total_energy_difference = total_energy + thermostat_work - total_energy.iloc[0]
  else:
    total_energy_difference = total_energy - total_energy.iloc[0]

  plt.figure()
  plt.plot(time_elapsed, total_energy_difference, label="Total energy error")

  diff = total_energy_difference[(time_elapsed > 0.7e-11) & (time_elapsed < 0.71e-11)]
  print(diff)

  plt.xlabel("time elapsed [s]")
  plt.ylabel("Energy [eV]")
  plt.title("Total energy error")
  # plt.xlim([2000000, 2100000])
  # plt.ylim([28, 31])
  plt.savefig(path + '/energy_difference.png')
  plt.show()

  plt.figure()
  plt.plot(time_elapsed, total_energy_show, label="Total energy")
  plt.xlabel("time elapsed [s]")
  plt.ylabel("Energy [eV]")
  plt.title("Total energy")
  plt.savefig(path + '/total_energy.png')
  plt.show()

  if thermostat:
    num_atoms = None
    with open(path + '/laamps/output_0.dump', 'r') as f:
      for line in f:
        if line.strip() == 'ITEM: NUMBER OF ATOMS':
          num_atoms = int(f.readline().strip())
          break

    if num_atoms is not None:
      mean_kinetic_energy = kinetic_energy / num_atoms
    else:
      mean_kinetic_energy = kinetic_energy  # Fallback if num_atoms not found

    T = 2 / 3 * mean_kinetic_energy / k_b

    plt.figure()
    plt.plot(time_elapsed, thermostat_epsilon, label="Thermostat epsilon")
    plt.xlabel("time elapsed [s]")
    plt.ylabel("thermostat epsilon")
    plt.title("Thermostat epsilon")
    plt.savefig(path + "/thermostat_epsilon.png")
    plt.show()

    plt.figure()
    plt.plot(time_elapsed, T, label="Temperature")
    plt.xlabel("time elapsed [s]")
    plt.ylabel("Temperature [K]")
    plt.title(f"Temperature")
    plt.ylim([0, 2400])
    plt.legend()
    plt.savefig(path + '/Temperature.png')
    plt.show()

    # print((T[15000:16000] - 2000).abs().max())
    # print((T[15000:16000] - 2000).abs().min())
    #
    # temp = T[15500:15800]
    # print(((temp - 2000).abs() < 15.0).sum())
    #
    # print((T - 1600).abs().min())


if __name__ == "__main__":
  # show_energy_plot_from_text("/home/washindeiru/studia/sem9/md/prog_check")
  # show_energy_plot_from_text("/home/washindeiru/studia/sem9/md/prog_check_v2/prog_check")
  import os
  # output_dir = "../../output"
  # output_dir = "/media/washindeiru/7E442D59442D1585/md"
  # newest_folder = max([os.path.join(output_dir, d) for d in os.listdir(output_dir)], key=os.path.getmtime)
  newest_folder = "/media/washindeiru/7E442D59442D1585/md/2025_05_29_more_dense_v2"
  newest_folder = "/media/washindeiru/7E442D59442D1585/md/2026_05_29_more_dense_v4"

  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/2026_05_29_more_dense_v4_with_lower_gravity"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/2026_05_29_more_dense_v2_with_lower_gravity"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/2026_05_29_more_dense_v4_with_1e-3_gravity"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/2025_05_29_more_dense_v4_from_begining_with_2e-3_gravity"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/2026_05_29_dense_v4_from_begining_with_lower_gravity_constant_temp"
  newest_folder = "/media/washindeiru/7E442D59442D1585/md/trash/results_cached/from_top_v2"
  newest_folder = "/media/washindeiru/7E442D59442D1585/md/plyn/from_top_v1_bigger_step"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/trash/results_cached/from_top_v2_07_06_2026_20_25_24"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/error_investigation/small_timestamp/z_edge"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/trash/results_cached/static_tube_results/2025_06_07_static_nanotube_test"
  newest_folder = "/media/washindeiru/7E442D59442D1585/md/trash/results_cached/from_top_v2"
  newest_folder = "/media/washindeiru/7E442D59442D1585/md/all/2026_05_29_more_dense_v4"
  newest_folder = "/media/washindeiru/7E442D59442D1585/md/plyn/from_top_v2_continued_local_v2"


  thermostat = True
  # newest_folder = "../../output/2026-04-14_12-12-07_exp"
  # compare_different_temps("../../output/2026-04-14_12-12-07_exp")
  show_energy_plot(newest_folder, thermostat, start=0, end=5e10)
  # compare_different_temps(newest_folder)
