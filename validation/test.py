import numpy as np
from numpy.ma.core import arccos
from collections import defaultdict
import csv
import matplotlib.pyplot as plt

def find_cosine(v1, v2):
    return np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))

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

if __name__ == '__main__':
    forces_data = read_forces_csv("../../output/2025-11-24_13-25-31/forces.csv")

    norm_1 = []
    norm_2 = []
    norm_3 = []

    cos_1_tab = []
    cos_2_tab = []
    cos_3_tab = []

    angle_1 = []
    angle_2 = []
    angle_3 = []

    for i in forces_data:
        v1 = forces_data[i][0]
        v2 = forces_data[i][1]
        v3 = forces_data[i][2]

        norm_1.append(np.linalg.norm(v1))
        norm_2.append(np.linalg.norm(v2))
        norm_3.append(np.linalg.norm(v3))

        cos1 = find_cosine(v1, v3)
        cos2 = find_cosine(v2, v3)
        cos3 = find_cosine(v2, v1)

        cos_1_tab.append(cos1)
        cos_2_tab.append(cos2)
        cos_3_tab.append(cos3)

        angle_1.append(arccos(cos1) * 180 / np.pi)
        angle_2.append(arccos(cos2) * 180 / np.pi)
        angle_3.append(arccos(cos3) * 180 / np.pi)
    # v1 = np.array([0.91276294686117521, 1.5809811135721659, 0])
    # v2 = np.array([0.91276294686117521, -1.5809811135721659, 0])
    # v3 = np.array([-1.8255675806274383, 0, 0])

    # v1 = np.array([-0.0163449510307766, 0, 0])
    # v2 = np.array([0.0081722888965992946, -0.014155082044696624, 0])
    # v3 = np.array([0.0081722888965992946, 0.014155082044696624, 0])

    plt.figure()
    plt.plot(norm_1, label="force 1")
    plt.plot(norm_2, label="force 2")
    plt.plot(norm_3, label="force 3")
    plt.title("długość wektora siły")
    plt.legend()
    plt.show()

    plt.figure()
    plt.plot(angle_1, label="angle 1")
    plt.title("Kąt pomiędzy siłą 1. i 2.")
    plt.show()

    plt.figure()
    plt.plot(angle_2, label="angle 2")
    plt.show()

    plt.figure()
    plt.plot(angle_3, label="angle 3")
    # plt.legend()
    plt.show()
    # plt.figure()
    # plt.plot(norm_1, label="Forces")
    # plt.plot(norm_2, label="Forces")

