#!/usr/bin/env python3
"""Particle-count scaling: fixed CPUs, vary num_particles via the random generator."""

from __future__ import annotations

import argparse
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path

SRC = Path("/mnt/7E442D59442D1585/md/paper/exp/particles-scaling/test")
BASE = Path("/mnt/7E442D59442D1585/md/paper/exp/particles-scaling")
BIN = Path("/home/washindeiru/studia/sem9/md/carbon_nanotube/target/release/carbon_nanotube")
REPEATS = [1]


def particle_counts() -> list[int]:
  """50, 100, 250, 500, then every 250 up to 7000."""
  # counts = [50, 100, 250, 500]
  counts = []
  n = 2250
  while n <= 7000:
    counts.append(n)
    n += 250
  return counts


def patch_parameters(params_path: Path, dest: Path) -> None:
  with params_path.open() as f:
    data = json.load(f)
  data["save_options"]["save_path"] = str(dest)
  data["save_options"]["keep_path"] = True
  with params_path.open("w") as f:
    json.dump(data, f, indent=2)
    f.write("\n")


def patch_generator(generator_path: Path, num_particles: int) -> None:
  with generator_path.open() as f:
    data = json.load(f)
  data["generator"]["config"]["num_particles"] = num_particles
  with generator_path.open("w") as f:
    json.dump(data, f, indent=2)
    f.write("\n")


def run_one(num_particles: int, repeat: int, dry_run: bool) -> None:
  dest = BASE / str(num_particles) / f"v{repeat}"
  print(f"=== N={num_particles} repeat=v{repeat} dest={dest} ===")

  if (dest / "energy.csv").is_file():
    print(f"skip existing {dest}")
    return

  if dry_run:
    print(f"dry-run: would run in {dest}")
    return

  dest.mkdir(parents=True, exist_ok=True)
  shutil.copy2(SRC / "parameters.json", dest / "parameters.json")
  shutil.copy2(SRC / "generator_config.json", dest / "generator_config.json")
  patch_parameters(dest / "parameters.json", dest)
  patch_generator(dest / "generator_config.json", num_particles)

  env = os.environ.copy()
  env["RUST_LOG"] = "info"
  subprocess.run(
    [str(BIN), "run", "--sg", str(dest)],
    env=env,
    check=True,
  )


def main() -> int:
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument(
    "--dry-run",
    action="store_true",
    help="Print planned runs without writing configs or starting the binary",
  )
  parser.add_argument(
    "--only",
    type=int,
    nargs="+",
    metavar="N",
    help="Run only these particle counts (default: full sweep)",
  )
  args = parser.parse_args()

  if not SRC.is_dir():
    print(f"template directory missing: {SRC}", file=sys.stderr)
    return 1
  if not (SRC / "parameters.json").is_file() or not (SRC / "generator_config.json").is_file():
    print(f"missing parameters.json or generator_config.json in {SRC}", file=sys.stderr)
    return 1
  if not args.dry_run and not BIN.is_file():
    print(f"binary missing: {BIN}", file=sys.stderr)
    return 1

  counts = args.only if args.only is not None else particle_counts()
  print(f"particle counts ({len(counts)}): {counts}")

  for n in counts:
    for j in REPEATS:
      run_one(n, j, dry_run=args.dry_run)

  print("=== all runs finished ===")
  return 0


if __name__ == "__main__":
  raise SystemExit(main())
