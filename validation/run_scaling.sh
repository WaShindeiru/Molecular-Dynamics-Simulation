#!/usr/bin/env bash
set -euo pipefail

SRC="/mnt/7E442D59442D1585/md/paper/exp/cores-scaling/result/good/big/test"
BASE="/mnt/7E442D59442D1585/md/paper/exp/cores-scaling/result/good/big"
BIN="/home/washindeiru/studia/sem9/md/carbon_nanotube/target/release/carbon_nanotube"

patch_save_path() {
  local params="$1"
  local dest="$2"
  python3 - "$params" "$dest" <<'PY'
import json, sys

path, dest = sys.argv[1], sys.argv[2]
with open(path) as f:
    data = json.load(f)
data["save_options"]["save_path"] = dest
data["save_options"]["keep_path"] = True
with open(path, "w") as f:
    json.dump(data, f, indent=2)
    f.write("\n")
PY
}

run_repeats() {
  local i="$1"
  local cpus="$2"
  for j in 1 2 3; do
    dest="$BASE/$i/v$j"
    echo "=== i=$i (taskset ${cpus}) repeat=v${j} dest=${dest} ==="
    if [[ -f "$dest/energy.csv" ]]; then
      echo "skip existing $dest"
      continue
    fi
    mkdir -p "$dest"
    cp "$SRC/parameters.json" "$SRC/generator_config.json" "$dest/"
    patch_save_path "$dest/parameters.json" "$dest"
    RUST_LOG=info taskset -c "$cpus" "$BIN" run --sg "$dest"
  done
}

run_repeats 11 "0-11"
# run_repeats 0 "0-0"

echo "=== all runs finished ==="
