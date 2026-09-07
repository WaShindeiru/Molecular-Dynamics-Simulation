import json

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import yaml
from pathlib import Path
from matplotlib.colors import to_rgb

import argparse
import os

from numpy.ma.core import size

TEMPERATURE_U = 11608.7

COLOR_SERIES = "#548b8d"

FONTSIZE_TITLE = 19
FONTSIZE_LABEL = 16
FONTSIZE_TICK = 13
FONTSIZE_LEGEND = 14

SCALING=1
FIGSIZE = (SCALING*6.5, SCALING*4.9)

# SCALING=0.832
# FIGSIZE = (SCALING*6.5, SCALING*4.9)
GRID_ALPHA = 0.3
LINEWIDTH = 2.2
LINEWIDTH_RAW = 1
LINEWIDTH_MEAN = 2.0
plt.rcParams["lines.linewidth"] = LINEWIDTH
DEFAULT_COLORS = plt.rcParams["axes.prop_cycle"].by_key()["color"]
MULTI_SERIES_COLORS = DEFAULT_COLORS[1:]


def series_mean_color():
  r, g, b = to_rgb(COLOR_SERIES)
  return (r * 0.65, g * 0.65, b * 0.65)


def new_figure(multi_series: bool = False) -> None:
  plt.figure(figsize=FIGSIZE)
  if multi_series:
    plt.gca().set_prop_cycle(color=MULTI_SERIES_COLORS)


def finish_plot(
  x_label: str,
  y_label: str,
  title: str,
  path: str,
  filename: str,
  legend: bool = True,
  legend_loc: str | None = None,
  ylim_top: float | None = None,
) -> None:
  plt.xlabel(x_label, fontsize=FONTSIZE_LABEL)
  plt.ylabel(y_label, fontsize=FONTSIZE_LABEL)
  plt.title(title, fontsize=FONTSIZE_TITLE)
  plt.grid(True, alpha=GRID_ALPHA)
  plt.tick_params(axis="both", labelsize=FONTSIZE_TICK)
  ax = plt.gca()
  ax.xaxis.get_offset_text().set_fontsize(FONTSIZE_LABEL)
  ax.yaxis.get_offset_text().set_fontsize(FONTSIZE_LABEL)
  if ylim_top is not None:
    plt.ylim(top=ylim_top)
  if legend:
    legend_kwargs = {"fontsize": FONTSIZE_LEGEND}
    if legend_loc is not None:
      legend_kwargs["loc"] = legend_loc
    plt.legend(**legend_kwargs)
  plt.tight_layout()
  plt.savefig(path + "/" + filename)
  plt.show()
  plt.close()


def load_time_step(path: str) -> float:
  parameters_path = Path(path) / "parameters.json"
  with open(parameters_path, encoding="utf-8") as f:
    parameters = json.load(f)
  return float(parameters["time_step"])


def plot_rolling_series(
  time_elapsed,
  values,
  x_label: str,
  y_label: str,
  title: str,
  path: str,
  filename: str,
  series_label: str,
  ylim_top: float | None = None,
  hline: float | None = None,
  hline_label: str | None = None,
) -> None:
  window = max(11, (len(values) // 400) | 1)
  roll = values.rolling(window=window, center=True, min_periods=1)
  smoothed = roll.mean()
  roll_min = roll.min()
  roll_max = roll.max()

  plt.figure(figsize=FIGSIZE)
  plt.fill_between(time_elapsed, roll_min, roll_max, color=COLOR_SERIES, alpha=0.22, label=series_label)
  plt.plot(time_elapsed, values, color=COLOR_SERIES, alpha=0.5, linewidth=LINEWIDTH_RAW)
  plt.plot(time_elapsed, smoothed, color=series_mean_color(), linewidth=LINEWIDTH_MEAN, label="Średnia ruchoma")
  if hline is not None:
    plt.axhline(hline, color="red", linestyle="--", linewidth=LINEWIDTH, label=hline_label)
  finish_plot(x_label, y_label, title, path, filename, ylim_top=ylim_top)


def show_energy_plot(path: str, use_time: bool = True, start: float | None = None, end: float | None = None, temp_ylim: float | None = None, cut: bool = False) -> None:
  energy_data = pd.read_csv(path + '/energy.csv', header=0)
  energy_data = energy_data[energy_data["iteration"] >= 1].reset_index(drop=True)

  iteration = energy_data["iteration"]
  if use_time:
    time_step = load_time_step(path)
    time_elapsed = iteration * time_step
    x_label = "czas [s]"
    if start is not None:
      energy_data = energy_data[time_elapsed >= start]
    if end is not None:
      energy_data = energy_data[time_elapsed <= end]
  else:
    x_label = "iteracja"
    if start is not None:
      energy_data = energy_data[iteration >= start]
    if end is not None:
      energy_data = energy_data[iteration <= end]
  energy_data = energy_data.reset_index(drop=True)

  has_thermostat = "thermostat_epsilon" in energy_data.columns
  has_nanotube_thermostat = "nanotube_thermostat_epsilon" in energy_data.columns

  iteration = energy_data["iteration"]
  if use_time:
    time_elapsed = iteration * time_step
  else:
    time_elapsed = iteration
  if cut and start is not None:
    time_elapsed = time_elapsed - start

  kinetic_energy_atom = energy_data["kinetic_energy_atom"]
  kinetic_energy_other = energy_data["kinetic_energy_other"]
  kinetic_energy = kinetic_energy_atom + kinetic_energy_other
  potential_energy = energy_data["potential_energy"]
  potential_gravity_energy = energy_data["potential_gravity_energy"]
  total_energy = energy_data["total_energy"]
  control_energy = energy_data["p_control_energy_total"]
  temperature = energy_data["temperature"] * TEMPERATURE_U
  if has_thermostat:
    thermostat_work_raw = energy_data["thermostat_work_total"]
    thermostat_epsilon = energy_data["thermostat_epsilon"]
    thermostat_work = thermostat_work_raw
    if cut and len(thermostat_work) > 0:
      thermostat_work = thermostat_work - thermostat_work.iloc[0]
  if has_nanotube_thermostat:
    nanotube_thermostat_epsilon = energy_data["nanotube_thermostat_epsilon"]
    nanotube_temperature = energy_data["nanotube_temperature"] * TEMPERATURE_U

  if has_thermostat:
    total_energy_show = thermostat_work + control_energy + kinetic_energy + potential_energy + potential_gravity_energy
  else:
    total_energy_show = control_energy + kinetic_energy + potential_energy + potential_gravity_energy

  new_figure(multi_series=True)
  plt.plot(time_elapsed, kinetic_energy, label="energia kinetyczna")
  plt.plot(time_elapsed, potential_energy, label="energia potencjalna")
  plt.plot(time_elapsed, total_energy_show, color=COLOR_SERIES, label="energia całkowita")
  # plt.plot(time_elapsed, control_energy, label="control energy")
  # plt.plot(time_elapsed, total_energy, label="all_energy")
  if has_thermostat:
    plt.plot(time_elapsed, thermostat_work, label="praca termostatu")
  # plt.xlim([130000, 140000])
  # plt.ylim(40000, 40100)
  # plt.xlim(500000, 630000)
  # plt.xlim([0, 3e-10])
  finish_plot(x_label, "Energia [eV]", "Energia", path, "energy.png")

  new_figure(multi_series=True)
  plt.plot(time_elapsed, kinetic_energy, label="energia kinetyczna")
  plt.plot(time_elapsed, potential_energy, label="energia potencjalna")
  plt.plot(time_elapsed, total_energy_show, color=COLOR_SERIES, label="energia całkowita")
  finish_plot(x_label, "Energia [eV]", "Energia", path, "energy_simple.png")

  new_figure(multi_series=True)
  plt.plot(time_elapsed, kinetic_energy, label="energia kinetyczna")
  plt.plot(time_elapsed, potential_energy, label="energia potencjalna")
  plt.plot(time_elapsed, total_energy_show, color=COLOR_SERIES, label="energia całkowita")
  plt.plot(time_elapsed, potential_gravity_energy, label="energia potencjalna\n grawitacji")
  finish_plot(x_label, "Energia [eV]", "Energia", path, "energy_simple_gravity.png", legend_loc="upper right")

  plt.figure(figsize=FIGSIZE)
  plt.plot(time_elapsed, kinetic_energy, color=COLOR_SERIES, label="energia kinetyczna")
  finish_plot(x_label, "Energia [eV]", "Energia kinetyczna", path, "kinetic_energy.png")
  plot_rolling_series(
    time_elapsed,
    kinetic_energy,
    x_label,
    "Energia [eV]",
    "Energia kinetyczna",
    path,
    "kinetic_energy_rolling.png",
    "energia kinetyczna",
  )

  plt.figure(figsize=FIGSIZE)
  plt.plot(time_elapsed, potential_energy, color=COLOR_SERIES, label="energia potencjalna")
  finish_plot(x_label, "Energia [eV]", "Energia potencjalna", path, "potential_energy.png")
  plot_rolling_series(
    time_elapsed,
    potential_energy,
    x_label,
    "Energia [eV]",
    "Energia potencjalna",
    path,
    "potential_energy_rolling.png",
    "energia potencjalna",
  )

  if has_thermostat:
    plt.figure(figsize=FIGSIZE)
    plt.plot(time_elapsed, thermostat_work, color=COLOR_SERIES, label="praca termostatu")
    finish_plot(x_label, "Energia [eV]", "Praca termostatu", path, "thermostat_work.png")
    plot_rolling_series(
      time_elapsed,
      thermostat_work,
      x_label,
      "Energia [eV]",
      "Praca termostatu",
      path,
      "thermostat_work_rolling.png",
      "praca termostatu",
    )

  plt.figure(figsize=FIGSIZE)
  plt.plot(time_elapsed, potential_gravity_energy, color=COLOR_SERIES, label="Energia potencjalna grawitacji")
  # plt.xlim([0, 3e-10])
  finish_plot(x_label, "Energia [eV]", "Energia potencjalna grawitacji", path, "gravitational_potential_energy.png")
  plot_rolling_series(
    time_elapsed,
    potential_gravity_energy,
    x_label,
    "Energia [eV]",
    "Energia potencjalna grawitacji",
    path,
    "gravitational_potential_energy_rolling.png",
    "energia potencjalna grawitacji",
  )

  if has_thermostat:
    total_energy_all = total_energy + thermostat_work_raw + control_energy
  else:
    total_energy_all = total_energy + control_energy

  # After --start/--end cropping, iloc[0] is the energy at the window start.
  energy_at_start = total_energy_all.iloc[0]
  total_energy_difference = total_energy_all - energy_at_start

  plt.figure(figsize=FIGSIZE)
  plt.plot(time_elapsed, total_energy_difference, color=COLOR_SERIES, label="Zmiana energii")

  max_error_label = total_energy_difference.diff().argmax()
  print(max_error_label)
  print(iteration[max_error_label])
  # diff = total_energy_difference[(time_elapsed > 0.7e-11) & (time_elapsed < 0.71e-11)]
  # print(diff)
  # plt.xlim([0, 2.3e-10])
  # plt.ylim([0, 60])
  # plt.xlim([0, 3e-10])
  finish_plot(x_label, "Energia [eV]", "Zmiana energii symulacji", path, "energy_difference.png")

  # Per-iteration absolute energy differences
  abs_energy_diff = np.abs(total_energy_all.diff().iloc[1:])
  cumulative_energy_change = abs_energy_diff.cumsum()

  # Create full series with first value = 0
  cumulative_energy_change_full = pd.Series(index=total_energy_all.index, dtype=float)
  cumulative_energy_change_full.iloc[0] = 0
  cumulative_energy_change_full.iloc[1:] = cumulative_energy_change.values

  plt.figure(figsize=FIGSIZE)
  plt.plot(time_elapsed, cumulative_energy_change_full, color=COLOR_SERIES, label="Błąd energii")
  finish_plot(x_label, "Energia [eV]", "Błąd energii symulacji", path, "energy_difference_momentum.png")

  plt.figure(figsize=FIGSIZE)
  plt.plot(
    time_elapsed,
    cumulative_energy_change_full / np.abs(energy_at_start),
    color=COLOR_SERIES,
    label="Względny błąd energii",
  )
  finish_plot(x_label, "Względny błąd energii", "Względny błąd energii symulacji", path, "energy_difference_momentum_relative.png")

  # plt.figure()
  # plt.plot(time_elapsed, total_energy_difference, label="Total energy error")

  # plt.xlabel(x_label)
  # plt.ylabel("Energy [eV]")
  # plt.title("Total energy error")
  # plt.xlim([0, 1e-10])
  # plt.ylim([-2, 20])
  # plt.savefig(path + '/energy_difference_small.png')
  # plt.show()

  plt.figure(figsize=FIGSIZE)
  plt.plot(time_elapsed, total_energy_show, color=COLOR_SERIES, label="Energia całkowita")
  # plt.xlim([0, 3e-10])
  finish_plot(x_label, "Energia [eV]", "Energia całkowita", path, "total_energy.png")

  plt.figure(figsize=FIGSIZE)
  plt.plot(time_elapsed, temperature, color=COLOR_SERIES, label="Temperature")
  # plt.axhline(1500, color="red", linestyle="--", label="Temperatura docelowa")
  # plt.xlim([0, 3e-10])
  finish_plot(x_label, "Temperatura [K]", "Temperatura", path, "Temperature.png", ylim_top=temp_ylim)

  plot_rolling_series(
    time_elapsed,
    temperature,
    x_label,
    "Temperatura [K]",
    "Temperatura",
    path,
    "Temperature_rolling.png",
    "Temperatura",
    ylim_top=temp_ylim,
  )

  if has_thermostat:
    plt.figure(figsize=FIGSIZE)
    # plt.plot(time_elapsed, thermostat_epsilon, color=COLOR_SERIES, label=r"$\xi$")
    plt.plot(time_elapsed, thermostat_epsilon, color=COLOR_SERIES, label="zmienna termostatu")
    # plt.xlim([0, 3e-10])
    # plt.ylim([-1, 1])
    finish_plot(x_label, r"zmienna termostatu $\xi$ [s$^{-1}$]", "Zmienna termostatu", path, "thermostat_epsilon.png")

    plot_rolling_series(
      time_elapsed,
      thermostat_epsilon,
      x_label,
      r"zmienna termostatu $\xi$ [s$^{-1}$]",
      "Zmienna termostatu",
      path,
      "thermostat_epsilon_rolling.png",
      r"$\xi$",
    )

    # print((T[15000:16000] - 2000).abs().max())
    # print((T[15000:16000] - 2000).abs().min())
    #
    # temp = T[15500:15800]
    # print(((temp - 2000).abs() < 15.0).sum())
    #
    # print((T - 1600).abs().min())

  if has_nanotube_thermostat:
    plt.figure(figsize=FIGSIZE)
    plt.plot(time_elapsed, nanotube_temperature, color=COLOR_SERIES, label="Nanotube temperature")
    finish_plot(
      x_label,
      "Temperatura [K]",
      "Temperatura nanorurki",
      path,
      "nanotube_temperature.png",
      ylim_top=temp_ylim,
    )

    plot_rolling_series(
      time_elapsed,
      nanotube_temperature,
      x_label,
      "Temperatura [K]",
      "Temperatura Nanorurki",
      path,
      "nanotube_temperature_rolling.png",
      "Temperatura",
      ylim_top=temp_ylim,
    )

    plt.figure(figsize=FIGSIZE)
    plt.plot(time_elapsed, nanotube_thermostat_epsilon, color=COLOR_SERIES, label="Nanotube thermostat epsilon")
    finish_plot(x_label, "thermostat epsilon", "Nanotube Thermostat epsilon", path, "nanotube_thermostat_epsilon.png")

    plot_rolling_series(
      time_elapsed,
      nanotube_thermostat_epsilon,
      x_label,
      "epsilon",
      "Epsilon termostatu nanorurki",
      path,
      "nanotube_thermostat_epsilon_rolling.png",
      "Epsilon",
    )

  evaluation = {"mean_temperature": float(temperature.mean())}
  if has_nanotube_thermostat:
    evaluation["mean_nanotube_temperature"] = float(nanotube_temperature.mean())

  with open(path + "/evaluation.yaml", "w", encoding="utf-8") as f:
    yaml.safe_dump(evaluation, f)


if __name__ == "__main__":

  parser = argparse.ArgumentParser(description="Plot energy data from a simulation output directory.")
  parser.add_argument("path", type=str, help="Path to the simulation output directory")
  parser.add_argument("--use-time", action=argparse.BooleanOptionalAction, default=True, help="Use time as x-axis (default: True)")
  parser.add_argument("--start", type=float, default=0, help="Start time [s] if --use-time, otherwise start iteration")
  parser.add_argument("--end", type=float, default=5e50, help="End time [s] if --use-time, otherwise end iteration")
  parser.add_argument("--temp-ylim", type=float, default=None, help="Upper y-axis limit for temperature plot")
  parser.add_argument("--cut", action="store_true", help="Shift time so --start is 0, and zero thermostat work at the first sample")
  args = parser.parse_args()

  # output_dir = "../../output"
  # output_dir = "/media/washindeiru/7E442D59442D1585/md"
  # newest_folder = max([os.path.join(output_dir, d) for d in os.listdir(output_dir)], key=os.path.getmtime)

  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/timestamp_investigation/trash/trash_v1_continued"
  # newest_folder = "/media/washindeiru/7E442D59442D1585/md/timestamp_investigation/triangle/e-17"

  # newest_folder = "../../output/2026-04-14_12-12-07_exp"
  # compare_different_temps("../../output/2026-04-14_12-12-07_exp")

  show_energy_plot(args.path, use_time=args.use_time, start=args.start, end=args.end, temp_ylim=args.temp_ylim, cut=args.cut)
  # compare_different_temps(args.path)
