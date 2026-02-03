import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def show_energy_plot(path: str, thermostat: bool) -> None:
  energy_data = pd.read_csv(path + '/energy.csv', header=0)
  iteration = energy_data.iloc[:, 0]
  kinetic_energy = energy_data.iloc[:, 1]
  potential_energy = energy_data.iloc[:, 2]
  total_energy = energy_data.iloc[:, 3]
  if thermostat:
    thermostat_work = energy_data.iloc[:, 4]
    thermostat_epsilon = energy_data.iloc[:, 5]

  potential_energies_per_atom = pd.read_csv(path + '/potential_energies.csv', header=None)
  potential_energies_per_atom = potential_energies_per_atom[potential_energies_per_atom.iloc[:, 0] != 1]
  potential_energy_summed = potential_energies_per_atom.groupby(potential_energies_per_atom.iloc[:, 0])[potential_energies_per_atom.columns[2]].sum()

  potential_energy_difference = potential_energy.values - potential_energy_summed.values

  plt.plot(iteration, kinetic_energy, label="kinetic energy")
  plt.plot(iteration, potential_energy, label="potential energy")
  plt.plot(iteration, total_energy, label="total energy")
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
  plt.plot(iteration, potential_energy_summed, label="potential energy summed")
  plt.xlabel("iteration")
  plt.ylabel("potential energy")
  plt.title("potential energy summed")
  plt.legend()
  plt.show()

  # Plot potential energy difference
  plt.figure()
  plt.plot(iteration, potential_energy_difference, label="Potential energy difference")
  plt.xlabel("iteration")
  plt.ylabel("Energy difference [eV]")
  plt.title("Potential energy difference (energy.csv vs summed per-atom)")
  plt.legend()
  plt.savefig(path + '/potential_energy_difference.png')
  plt.show()

  # Plot potential energy of the first particle in every iteration
  first_particle_id = 0
  first_particle_energy = potential_energies_per_atom[potential_energies_per_atom.iloc[:, 1] == first_particle_id]

  positions = pd.read_csv(path + '/positions.csv', header=None, names=['iteration', 'particle_id', 'x', 'y', 'z'])
  positions = positions[positions.iloc[:, 0] != 1]
  first_particle_positions = positions[positions['particle_id'] == first_particle_id]
  initial_x = first_particle_positions['x'].iloc[0]
  x_displacement = first_particle_positions['x'] - initial_x

  plt.figure()
  plt.plot(x_displacement, first_particle_energy.iloc[:, 2], label=f"Potential energy of particle {first_particle_id}")
  plt.xlabel("iteration")
  plt.ylabel("Potential energy [eV]")
  plt.title(f"Potential energy of particle {first_particle_id}")
  plt.legend()
  plt.xlim([0, 4])
  plt.ylim([-1, 4])
  plt.savefig(path + '/first_particle_potential_energy.png')
  plt.show()

  # Read positions and calculate x displacement of first particle


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
    plt.figure()
    plt.plot(iteration, thermostat_epsilon, label="Thermostat epsilon")
    plt.xlabel("iteration")
    plt.ylabel("thermostat epsilon")
    plt.title("Thermostat epsilon")
    plt.savefig(path + "/thermostat_epsilon.png")
    plt.show()


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
  thermostat = False
  show_energy_plot(newest_folder, thermostat)
