import json
import numpy as np
from enum import Enum
from dataclasses import dataclass, field


class AtomType(str, Enum):
  C_nanotube = "C_nanotube"
  Fe_nanotube = "Fe_nanotube"


class ValueUnits(str, Enum):
  si = "si"
  unitless = "unitless"


@dataclass
class Position:
  x: float
  y: float
  z: float

  @staticmethod
  def from_numpy(v: np.ndarray) -> "Position":
    return Position(x=float(v[0]), y=float(v[1]), z=float(v[2]))

  def to_numpy(self) -> np.ndarray:
    return np.array([self.x, self.y, self.z])


@dataclass
class Particle:
  atom_type: AtomType
  position: Position


@dataclass
class Simulation:
  particles: list[Particle]
  value_units: ValueUnits = ValueUnits.unitless
  offset: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
  vel_mean: list[float] = field(default_factory=lambda: [0.0, 0.0, 0.0])
  vel_std_dev: list[float] = field(default_factory=lambda: [1e-3, 1e-3, 1e-3])

  def save_json(self, filepath: str, indent: int = 2) -> None:
    with open(filepath, "w") as f:
      f.write(self.to_json(indent=indent))

  def to_json(self, indent: int = 2) -> str:
    data = {
      "value_units": self.value_units.value,
      "generator": {
        "type": "nanotube",
        "config": {
          "particles": [
            {
              "position": {"x": p.position.x, "y": p.position.y, "z": p.position.z},
              "particle_type": p.atom_type.value,
            }
            for p in self.particles
          ],
          "offset": {"x": self.offset[0], "y": self.offset[1], "z": self.offset[2]},
          "vel_mean": {"x": self.vel_mean[0], "y": self.vel_mean[1], "z": self.vel_mean[2]},
          "vel_std_dev": {"x": self.vel_std_dev[0], "y": self.vel_std_dev[1], "z": self.vel_std_dev[2]},
        },
      },
    }
    return json.dumps(data, indent=indent)


_XYZ_ATOM_MAP: dict[str, AtomType] = {
  "C": AtomType.C_nanotube,
  "Fe": AtomType.Fe_nanotube,
}


def parse_xyz(filepath: str) -> Simulation:
  with open(filepath) as f:
    lines = f.readlines()

  n_atoms = int(lines[0].strip())

  particles = []
  for line in lines[1:1 + n_atoms]:
    parts = line.split()
    atom_type = _XYZ_ATOM_MAP.get(parts[0], AtomType.C_nanotube)
    position = Position(x=float(parts[1]), y=float(parts[2]), z=float(parts[3]))
    particles.append(Particle(atom_type=atom_type, position=position))

  return Simulation(particles=particles)


if __name__ == "__main__":
  sim = parse_xyz("/home/washindeiru/studia/sem9/md/carbon_tubes/minimum_tube_v2_27_05_26/Nanotube_CC.txt")
  print(f"number of atoms: {len(sim.particles)}")
  sim.save_json("/home/washindeiru/studia/sem9/md/carbon_tubes/minimum_tube_v2_27_05_26/Nanotube_CC.json")
  print("saved to nanotube.json")
