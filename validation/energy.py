import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def show_energy_plot(path: str) -> None:
    energy_data = pd.read_csv(path + '/energy.csv', header=1)
    iteration = energy_data.iloc[:, 0]
    kinetic_energy = energy_data.iloc[:, 1]
    potential_energy = energy_data.iloc[:, 2]
    total_energy = energy_data.iloc[:, 3]

    plt.plot(iteration, kinetic_energy, label="kinetic energy")
    plt.plot(iteration, potential_energy, label="potential energy")
    plt.plot(iteration, total_energy, label="total energy")

    plt.xlabel("iteration")
    plt.ylabel("Energy [eV]")
    # plt.xlim([0, 1000])
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
    # plt.xlim([0, 1000])
    plt.savefig(path + '/energy_difference.png')
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
    show_energy_plot("../../output/2026-01-18_23-58-52")