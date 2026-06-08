import re
import shutil
from pathlib import Path

_DUMP_PATTERN = re.compile(r"^output_(\d+)\.dump$")


def _dump_iteration(path: Path) -> int | None:
  match = _DUMP_PATTERN.match(path.name)
  if match is None:
    return None
  return int(match.group(1))


def merge(first_path: str, second_path: str, iter_num: int, output_path: str) -> None:
  """Merge LAMMPS dump files from two directories into ``output_path``.

  - From ``first_path``: copies ``output_{i}.dump`` for all ``i < iter_num``.
  - From ``second_path``: copies every ``output_{j}.dump`` as
    ``output_{j + iter_num}.dump`` in the output directory.

  Iteration numbering starts at 0.
  """
  first_dir = Path(first_path)
  second_dir = Path(second_path)
  output_dir = Path(output_path)

  if not first_dir.is_dir():
    raise FileNotFoundError(f"Directory does not exist: {first_path}")
  if not second_dir.is_dir():
    raise FileNotFoundError(f"Directory does not exist: {second_path}")

  output_dir.mkdir(parents=True, exist_ok=True)

  for dump_file in first_dir.iterdir():
    if not dump_file.is_file():
      continue
    iteration = _dump_iteration(dump_file)
    if iteration is None or iteration >= iter_num:
      continue
    shutil.copy2(dump_file, output_dir / dump_file.name)

  for dump_file in second_dir.iterdir():
    if not dump_file.is_file():
      continue
    iteration = _dump_iteration(dump_file)
    if iteration is None:
      continue
    shifted_name = f"output_{iteration + iter_num}.dump"
    shutil.copy2(dump_file, output_dir / shifted_name)


if __name__ == "__main__":
  merge("/media/washindeiru/7E442D59442D1585/md/trash/results_cached/from_top_v2",
        "/media/washindeiru/7E442D59442D1585/md/plyn/from_top_v2_continued_local_v2/laamps",
        1847771,
        "/media/washindeiru/7E442D59442D1585/md/plyn/from_top_v2_merged/laamps")