#!/usr/bin/env python3
import argparse
import re
from pathlib import Path

import pandas as pd

TIME_PATTERN = re.compile(
  r"^Simulation Time:\s*(?P<seconds>[0-9.]+)s\s+seconds\s*$",
  re.MULTILINE,
)
CORE_DIR_PATTERN = re.compile(r"^\d+$")
VERSION_DIR_PATTERN = re.compile(r"^v(\d+)$")


def parse_simulation_time(info_path: Path) -> float:
  text = info_path.read_text(encoding="utf-8")
  match = TIME_PATTERN.search(text)
  if match is None:
    raise ValueError(f"Could not parse simulation time in {info_path}")
  return float(match.group("seconds"))


def load_core_times(directory: str | Path) -> pd.DataFrame:
  root = Path(directory)
  if not root.is_dir():
    raise FileNotFoundError(f"Directory does not exist: {directory}")

  records: dict[int, dict[str, float]] = {}
  versions: set[str] = set()

  for core_dir in sorted(root.iterdir(), key=lambda p: (not p.name.isdigit(), int(p.name) if p.name.isdigit() else p.name)):
    if not core_dir.is_dir() or CORE_DIR_PATTERN.match(core_dir.name) is None:
      continue

    core = int(core_dir.name)
    for version_dir in core_dir.iterdir():
      if not version_dir.is_dir():
        continue
      version_match = VERSION_DIR_PATTERN.match(version_dir.name)
      if version_match is None:
        continue

      info_path = version_dir / "info.txt"
      if not info_path.is_file():
        continue

      version = version_dir.name
      versions.add(version)
      records.setdefault(core, {})[version] = parse_simulation_time(info_path)

  if not records:
    raise FileNotFoundError(f"No core/version info.txt files found in: {directory}")

  version_columns = sorted(versions, key=lambda name: int(VERSION_DIR_PATTERN.match(name).group(1)))
  rows = []
  for core in sorted(records):
    row = {"cores": core}
    times = []
    for version in version_columns:
      time = records[core].get(version)
      row[version] = time
      if time is not None:
        times.append(time)
    row["mean"] = sum(times) / len(times) if times else float("nan")
    rows.append(row)

  return pd.DataFrame(rows, columns=["cores", *version_columns, "mean"])


def main() -> None:
  parser = argparse.ArgumentParser(
    description="Build a cores × experiment-time dataframe from info.txt files.",
  )
  parser.add_argument("directory", help="Directory containing 0..11 subdirectories")
  args = parser.parse_args()

  dataframe = load_core_times(args.directory)
  pd.set_option("display.max_columns", None)
  pd.set_option("display.width", None)
  print(dataframe.to_string(index=False))

  dataframe.to_parquet("./cpu.parquet")


if __name__ == "__main__":
  main()
