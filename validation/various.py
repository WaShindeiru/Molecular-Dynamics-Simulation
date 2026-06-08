import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path

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

def get_distribution(data, idx_01, idx_02, name):
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
  plt.savefig(name)
  plt.show()

  return df_merged['distance'].std()

def compare_different_temps(path_: str):
  idx_01 = 10
  idx_02 = 40
  result = read_atoms_from_output_dumps(path_, idx_01, idx_02)
  temp = result[(result['iteration'] >= 19756) & (result['iteration'] < 29756)]

  deviation_1 = get_distribution(temp, idx_01, idx_02, path_ + "/first")
  print(deviation_1)

  temp_2 = result[(result['iteration'] >= 70000) & (result['iteration'] < 80000)]

  deviation_2 = get_distribution(temp_2, idx_01, idx_02, path_ + "/second")
  print(deviation_2)

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

  iteration = energy_data.iloc[1:, 0]
  potential_energy = energy_data.iloc[1:, 2]  # 3rd column (index 2)
  kinetic_energy = energy_data.iloc[1:, 3]  # 4th column (index 3)
  total_energy = energy_data.iloc[1:, 4]  # 5th column (index 4)

  plt.plot(iteration, kinetic_energy, label="kinetic energy")
  plt.plot(iteration, potential_energy, label="potential energy")
  plt.plot(iteration, total_energy, label="total energy")

  plt.xlabel("iteration")
  plt.ylabel("Energy [eV]")
  plt.title("Energy plot")

  plt.legend()
  plt.show()

  total_energy_difference = total_energy - total_energy.iloc[1]
  plt.figure()
  plt.plot(iteration, total_energy_difference, label="total energy difference")

  plt.xlabel("iteration")
  plt.ylabel("Energy [eV]")
  plt.title("Total energy difference")
  plt.show()