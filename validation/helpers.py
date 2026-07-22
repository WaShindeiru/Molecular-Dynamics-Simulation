import json
import sys


def transform_particles(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    for particle in data["particles"]:
        if (
            particle.get("particle_type") == "CustomVelocityAtom"
            and particle.get("velocity_manager_id") == 1
        ):
            particle["particle_type"] = "VelocityControlledParticle"
            del particle["velocity_manager_id"]
            particle["control_velocity_manager_id"] = 0
            particle["position"]["z"] += 4e-11

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def reset_velocity_manager_id(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    for particle in data["particles"]:
        if particle.get("velocity_manager_id") == 2:
            particle["velocity_manager_id"] = 0

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def move_custom_velocity_atoms(input_path, output_path, dz):
    with open(input_path) as f:
        data = json.load(f)

    for particle in data["particles"]:
        if particle.get("particle_type") == "VelocityControlledParticle":
            particle["position"]["z"] += dz

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def revert_high_velocity_controlled_particles(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    for particle in data["particles"]:
        if (
            particle.get("particle_type") == "VelocityControlledParticle"
            and particle["position"]["z"] > 4.5e-9
        ):
            assert particle.get("atom_type") == "C_nanotube", (
                f"expected atom_type 'C_nanotube', got {particle.get('atom_type')!r}"
            )
            particle["particle_type"] = "CustomVelocityAtom"
            del particle["control_velocity_manager_id"]
            particle["velocity_manager_id"] = 1
            particle["atom_type"] = "C_nanotube_static"

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


if __name__ == "__main__":
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    revert_high_velocity_controlled_particles(input_path, output_path)
    # transform_particles(input_path, output_path)
    # move_custom_velocity_atoms(input_path, output_path, 1e-11)
