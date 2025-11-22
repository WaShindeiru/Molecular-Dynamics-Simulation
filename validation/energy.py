import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def show_energy_plot(path: str) -> None:
    energy_data = pd.read_csv(path + '/energy.csv', header=None)
    iteration = energy_data.iloc[:, 0]
    kinetic_energy = energy_data.iloc[:, 1]
    potential_energy = energy_data.iloc[:, 2]
    total_energy = energy_data.iloc[:, 3]

    plt.plot(iteration, kinetic_energy, label="kinetic energy")
    plt.plot(iteration, potential_energy, label="potential energy")
    plt.plot(iteration, total_energy, label="total energy")

    plt.xlabel("iteration")
    plt.ylabel("Energy [eV]")
    plt.xlim([0, 1000])
    plt.title("Energy plot")

    plt.legend()
    plt.savefig(path + '/energy.png')
    plt.show()

    total_energy_difference = total_energy - total_energy[0]
    plt.figure()
    plt.plot(iteration, total_energy_difference, label="total energy difference")

    plt.xlabel("iteration")
    plt.ylabel("Energy [eV]")
    plt.title("Total energy difference")
    plt.xlim([0, 1000])
    plt.show()


if __name__ == "__main__":
    show_energy_plot("../../output/2025-11-22_00-32-45_simple_vel_verlet_e_16_single_atom")