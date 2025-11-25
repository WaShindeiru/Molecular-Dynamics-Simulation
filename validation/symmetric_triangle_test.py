import csv
import numpy as np
from collections import defaultdict

from matplotlib import pyplot as plt


def read_forces_csv(filename):
    forces_data = defaultdict(dict)

    with open(filename, 'r') as file:
        csv_reader = csv.reader(file)

        # Skip header if exists, otherwise remove this line
        # next(csv_reader, None)

        for row in csv_reader:
            if len(row) < 5:
                continue  # Skip incomplete rows

            iteration = int(row[0])
            atom_id = int(row[1])
            force_x = float(row[2])
            force_y = float(row[3])
            force_z = float(row[4])

            forces_data[iteration][atom_id] = np.array([force_x, force_y, force_z])

    return dict(forces_data)


if __name__ == "__main__":
    forces_data = read_forces_csv("../../output/symmetric_triangle_test/forces.csv")

    x_difference = []
    y_difference = []
    z_aha = []

    third_atom_y_error = []
    third_atom_z_error = []
    third_atom_x_force = []

    for i in range(1,100):
        x_diff = abs(forces_data[i][0][0] - forces_data[i][1][0])
        y_diff = abs(-forces_data[i][0][1] - forces_data[i][1][1])
        z_temp = abs(forces_data[i][0][2]) + abs( forces_data[i][1][2])

        x_difference.append(x_diff)
        y_difference.append(y_diff)
        z_aha.append(z_temp)

        third_atom_y_error.append(forces_data[i][2][1])
        third_atom_z_error.append(forces_data[i][2][2])
        third_atom_x_force.append(forces_data[i][2][0])

    index = np.linspace(0, len(x_difference), len(x_difference))
    plt.plot(index, x_difference)
    plt.show()
    plt.plot(index, y_difference)
    plt.show()
    plt.plot(index, z_aha)
    plt.show()

    plt.plot(index, third_atom_y_error)
    plt.show()
    plt.plot(index, third_atom_z_error)
    plt.show()
    plt.plot(index, third_atom_x_force)
    plt.show()