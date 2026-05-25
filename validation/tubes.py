import numpy as np

def parse_xyz(filepath):
    with open(filepath) as f:
        lines = f.readlines()

    n_atoms = int(lines[0].strip())

    atoms = []
    coords = []

    for line in lines[1:1 + n_atoms]:
        parts = line.split()
        atoms.append(parts[0])
        coords.append([float(parts[1]), float(parts[2]), float(parts[3])])

    return n_atoms, atoms, np.array(coords)


if __name__ == "__main__":
  n_atoms, atoms, coords = parse_xyz("/home/washindeiru/studia/sem9/md/carbon_tubes/Nanotube_CC.txt")
  print(f"number of atoms: {n_atoms}")
  print(atoms)   # ['C', 'C', 'C', ...]
  print(coords.shape)  # shape (120, 3), x/y/z for each atom
  print(coords)