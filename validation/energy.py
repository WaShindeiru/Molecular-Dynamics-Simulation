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
    plt.title("Energy plot")

    plt.legend()
    plt.show()


if __name__ == "__main__":
    show_energy_plot("../output/2025-11-14_20-41-00_semi_implitic_euler_e_16")