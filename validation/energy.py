import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

k_b = 8.617333e-5

def read_atoms_from_output_dumps(directory: str, atom_id_1: int, atom_id_2: int) -> pd.DataFrame:
  """Read positions of two selected atoms from all output*.dump files in a directory."""
  directory_path = Path(directory)
  if not directory_path.is_dir():
    raise FileNotFoundError(f"Directory does not exist: {directory}")

  dump_files = sorted(
    [p for p in directory_path.iterdir() if p.is_file() and p.name.startswith("output") and p.name.endswith("dump")]
  )

  if not dump_files:
    raise FileNotFoundError(f"No files matching output*.dump in: {directory}")

  selected_ids = {int(atom_id_1), int(atom_id_2)}
  records = []

  for dump_file in dump_files:
    with dump_file.open("r") as f:
      lines = f.readlines()

    i = 0
    current_iteration = None

    while i < len(lines):
      line = lines[i].strip()

      if line == "ITEM: TIMESTEP":
        i += 1
        if i < len(lines):
          current_iteration = int(lines[i].strip())
        i += 1
        continue

      if line.startswith("ITEM: ATOMS"):
        headers = line.split()[2:]
        required_headers = ("id", "x", "y", "z")
        if not all(h in headers for h in required_headers):
          raise ValueError(
            f"Missing one of {required_headers} in ATOMS header of file: {dump_file}"
          )

        id_idx = headers.index("id")
        x_idx = headers.index("x")
        y_idx = headers.index("y")
        z_idx = headers.index("z")
        max_idx = max(id_idx, x_idx, y_idx, z_idx)

        i += 1
        while i < len(lines) and not lines[i].startswith("ITEM:"):
          parts = lines[i].split()
          if len(parts) > max_idx:
            atom_id = int(float(parts[id_idx]))
            if atom_id in selected_ids:
              records.append(
                {
                  "iteration": current_iteration,
                  "atom_id": atom_id,
                  "x": float(parts[x_idx]),
                  "y": float(parts[y_idx]),
                  "z": float(parts[z_idx]),
                }
              )
          i += 1
        continue

      i += 1

  if not records:
    return pd.DataFrame(columns=["iteration", "atom_id", "x", "y", "z"])

  result = pd.DataFrame(records, columns=["iteration", "atom_id", "x", "y", "z"])
  return result.sort_values(["iteration", "atom_id"]).reset_index(drop=True)

def get_distribution(data, idx_01, idx_02):
  df_pivot = data.pivot_table(index='iteration', columns='atom_id')
  df_merged = pd.DataFrame({
    'iteration': df_pivot.index,
    'x_20': df_pivot['x'][idx_01].values,
    'y_20': df_pivot['y'][idx_01].values,
    'z_20': df_pivot['z'][idx_01].values,
    'x_30': df_pivot['x'][idx_02].values,
    'y_30': df_pivot['y'][idx_02].values,
    'z_30': df_pivot['z'][idx_02].values,
  })

  df_merged['distance'] = np.sqrt(
    (df_merged['x_30'] - df_merged['x_20']) ** 2 +
    (df_merged['y_30'] - df_merged['y_20']) ** 2 +
    (df_merged['z_30'] - df_merged['z_20']) ** 2
  )

  print(df_merged.shape)

  plt.figure()
  plt.plot(df_merged['iteration'], df_merged['distance'])
  plt.show()

  return df_merged['distance'].std()

def compare_different_temps(path_: str):
  idx_01 = 16
  idx_02 = 27
  result = read_atoms_from_output_dumps(path_, idx_01, idx_02)
  temp = result[(result['iteration'] >= 1218) & (result['iteration'] < 31218)]

  deviation_1 = get_distribution(temp, idx_01, idx_02)
  print(deviation_1)

  temp_2 = result[(result['iteration'] >= 50000) & (result['iteration'] < 80000)]

  deviation_2 = get_distribution(temp_2, idx_01, idx_02)
  print(deviation_2)


def show_energy_plot(path: str, thermostat: bool) -> None:
  energy_data = pd.read_csv(path + '/energy.csv', header=0)
  iteration = energy_data.iloc[:, 0]
  kinetic_energy = energy_data.iloc[:, 1]
  potential_energy = energy_data.iloc[:, 2]
  potential_gravity_energy = energy_data.iloc[:, 3]
  total_energy = energy_data.iloc[:, 4]
  if thermostat:
    thermostat_work = energy_data.iloc[:, 5]
    thermostat_epsilon = energy_data.iloc[:, 6]

  total_energy_show = kinetic_energy + potential_energy

  plt.plot(iteration, kinetic_energy, label="kinetic energy")
  plt.plot(iteration, potential_energy, label="potential energy")
  plt.plot(iteration, total_energy_show, label="total energy")
  if thermostat:
    plt.plot(iteration, thermostat_work, label="thermostat work")
  plt.xlabel("iteration")
  plt.ylabel("Energy [eV]")
  # plt.xlim([0, 1000])
  plt.title("Energy plot")
  plt.legend()
  plt.savefig(path + '/energy.png')
  plt.show()

  plt.figure()
  plt.plot(iteration, potential_gravity_energy, label="gravitational pot energy")
  plt.xlabel("iteration")
  plt.ylabel("Energy [eV]")
  plt.title("Gravitational Potential Energy")
  plt.legend()
  plt.savefig(path + '/gravitational_potential_energy.png')
  plt.show()

  if thermostat:
    total_energy_difference = total_energy + thermostat_work - total_energy[0]
  else:
    total_energy_difference = total_energy - total_energy[0]

  plt.figure()
  plt.plot(iteration, total_energy_difference, label="Total energy error")

  plt.xlabel("iteration")
  plt.ylabel("Energy [eV]")
  plt.title("Total energy error")
  # plt.xlim([0, 1000])
  plt.savefig(path + '/energy_difference.png')
  plt.show()

  if thermostat:
    num_atoms = None
    with open(path + '/output_0.dump', 'r') as f:
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
    plt.plot(iteration, thermostat_epsilon, label="Thermostat epsilon")
    plt.xlabel("iteration")
    plt.ylabel("thermostat epsilon")
    plt.title("Thermostat epsilon")
    plt.savefig(path + "/thermostat_epsilon.png")
    plt.show()

    plt.figure()
    plt.plot(iteration, T, label="Temperature")
    plt.xlabel("iteration")
    plt.ylabel("Temperature [K]")
    plt.title(f"Temperature")
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


def show_energy_plot_from_text(path: str) -> None:
  # Read fixed-width file with 5 columns, each 20 characters wide
  energy_data = pd.read_fwf(
    path + '/en.dat',
    widths=[21, 21, 21, 21, 21],
    header=None
  )

  # Clean and convert to numeric, removing any non-printable characters
  # for col in energy_data.columns:
  #     energy_data[col] = pd.to_numeric(energy_data[col].astype(str).str.strip(), errors='coerce')

  iteration = energy_data.iloc[:, 0]
  potential_energy = energy_data.iloc[:, 2]  # 3rd column (index 2)
  kinetic_energy = energy_data.iloc[:, 3]  # 4th column (index 3)
  total_energy = energy_data.iloc[:, 4]  # 5th column (index 4)

  plt.plot(iteration, kinetic_energy, label="kinetic energy")
  plt.plot(iteration, potential_energy, label="potential energy")
  plt.plot(iteration, total_energy, label="total energy")

  plt.xlabel("iteration")
  plt.ylabel("Energy [eV]")
  plt.title("Energy plot")

  plt.legend()
  plt.show()

  total_energy_difference = total_energy - total_energy.iloc[0]
  plt.figure()
  plt.plot(iteration, total_energy_difference, label="total energy difference")

  plt.xlabel("iteration")
  plt.ylabel("Energy [eV]")
  plt.title("Total energy difference")
  plt.show()

if __name__ == "__main__":
  # show_energy_plot_from_text("/home/washindeiru/studia/sem9/md/prog_check")
  # show_energy_plot_from_text("/home/washindeiru/studia/sem9/md/prog_check_v2/prog_check")
  import os
  output_dir = "../../output"
  newest_folder = max([os.path.join(output_dir, d) for d in os.listdir(output_dir)], key=os.path.getmtime)
  thermostat = True
  # compare_different_temps("../../output/2026-03-31_00-08-46")
  # show_energy_plot(newest_folder, thermostat)
  compare_different_temps(newest_folder)
