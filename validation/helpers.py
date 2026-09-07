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


def remove_nanotube_particles(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    data["particles"] = [
        particle
        for particle in data["particles"]
        if particle.get("atom_type") not in ("C_nanotube_static", "C_nanotube")
    ]

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def keep_atom_particles(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    data["particles"] = [
        particle
        for particle in data["particles"]
        if particle.get("particle_type") == "Atom"
    ]

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def convert_velocity_controlled_into_atom(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    for particle in data["particles"]:
        if particle["particle_type"] == "VelocityControlledParticle":
            assert("control_velocity_manager_id" in particle)
            particle["particle_type"] = "Atom"
            del particle["control_velocity_manager_id"]

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def filter_particles_in_box(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    data["particles"] = [
        particle
        for particle in data["particles"]
        if (
            0 <= particle["position"]["x"] <= 100e-10
            and 0 <= particle["position"]["y"] <= 100e-10
            and 0 <= particle["position"]["z"] <= 14e-10
        )
    ]


    for i, particle in enumerate(data["particles"]):
        particle["id"] = i

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def move_particles_by_type(input_path, output_path, particle_type, dx=0.0, dy=0.0, dz=0.0):
    with open(input_path) as f:
        data = json.load(f)

    for particle in data["particles"]:
        if particle.get("particle_type") == particle_type:
            particle["position"]["x"] += dx
            particle["position"]["y"] += dy
            particle["position"]["z"] += dz

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def translate_particles(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    for particle in data["particles"]:
        particle["position"]["x"] += 9.2e-10
        particle["position"]["y"] += 9.55e-10
        particle["position"]["z"] += 2.47e-9

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def freeze_low_z_particles(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    for particle in data["particles"]:
        if particle["position"]["z"] > 20e-10:
            particle["atom_type"] = "C_nanotube_static"
            particle["particle_type"] = "CustomVelocityAtom"
            particle["velocity_manager_id"] = 1
            particle.pop("control_velocity_manager_id", None)

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


CARBON_ATOM_TYPES = ("C", "C_nanotube", "C_nanotube_static")


def reindex_particles(input_path, output_path):
    with open(input_path) as f:
        data = json.load(f)

    carbon = 0
    iron = 0
    for i, particle in enumerate(data["particles"]):
        particle["id"] = i
        atom_type = particle.get("atom_type")
        if atom_type in CARBON_ATOM_TYPES:
            carbon += 1
        elif atom_type == "Fe":
            iron += 1
        else:
            raise ValueError(f"unknown atom_type {atom_type!r}")

    data["num_of_atoms"] = len(data["particles"])
    data["num_of_carbon_atoms"] = carbon
    data["num_of_iron_atoms"] = iron

    with open(output_path, "w") as f:
        json.dump(data, f, indent=2)


def inspect_particle_bounds(input_path):
    with open(input_path) as f:
        data = json.load(f)

    particles = data["particles"]
    if not particles:
        raise ValueError(f"no particles in {input_path}")

    xs = [particle["position"]["x"] for particle in particles]
    ys = [particle["position"]["y"] for particle in particles]
    zs = [particle["position"]["z"] for particle in particles]

    print(f"x: min = {min(xs)}, max = {max(xs)}")
    print(f"y: min = {min(ys)}, max = {max(ys)}")
    print(f"z: min = {min(zs)}, max = {max(zs)}")


if __name__ == "__main__":
    input_path = sys.argv[1]
    output_path = sys.argv[2]
    move_particles_by_type(
        input_path,
        output_path,
        "Atom",
        dx=0.0,
        dy=0.0,
        dz=1e-9,
    )
    # reindex_particles(input_path, output_path)
    # reindex_particles(input_path, output_path)
    # translate_particles(input_path, output_path)
    # freeze_low_z_particles(input_path, output_path)
    # remove_nanotube_particles(input_path, output_path)
    # revert_high_velocity_controlled_particles(input_path, output_path)
    # transform_particles(input_path, output_path)
    # move_custom_velocity_atoms(input_path, output_path, 1e-11)
