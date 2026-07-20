#!/bin/bash

# Runs RadiSeq concurrently, once per simulation_* folder in run_results/ that has a
# parameters.txt but does not yet have an output folder (as declared by that folder's own
# output_directory_path parameter, defaulting to ./output if unset). Folders that already
# have their output are skipped, so this can be re-run to pick up only new/incomplete
# simulations. Each run's stdout/stderr goes to run.log inside its own folder, since running
# several RadiSeq processes at once would otherwise interleave their terminal output. Waits
# for every launched run to finish before exiting.

script_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo_root="$(cd "$script_dir/.." && pwd)"
run_results_dir="$repo_root/run_results"
radiseq_bin="$repo_root/RadiSeq"

export RADISEQ_DATA_DIR="${RADISEQ_DATA_DIR:-$repo_root/radiSeqData}"

if [ ! -x "$radiseq_bin" ]; then
    echo "RadiSeq binary not found or not executable: $radiseq_bin"
    exit 1
fi

if [ ! -d "$run_results_dir" ]; then
    echo "run_results folder not found: $run_results_dir"
    exit 1
fi

pids=()
names=()

for sim_dir in "$run_results_dir"/*/; do
    sim_dir="${sim_dir%/}"
    sim_name="$(basename "$sim_dir")"

    [ -f "$sim_dir/parameters.txt" ] || continue

    output_rel=$(grep -E "^[[:space:]]*output_directory_path[[:space:]]*=" "$sim_dir/parameters.txt" \
        | head -n1 | sed -E 's/^[^=]*=[[:space:]]*//; s/[[:space:]]*#.*//; s/[[:space:]]*$//')
    output_rel="${output_rel:-./output}"
    output_dir="$sim_dir/${output_rel#./}"

    if [ -d "$output_dir" ]; then
        echo "Skipping $sim_name (output already exists)"
        continue
    fi

    if [ -f "$sim_dir/description.json" ]; then
        python3 -c "
import json
from datetime import datetime

path = '$sim_dir/description.json'
with open(path) as f:
    description = json.load(f)
description['time_simulation_run'] = datetime.now().isoformat()
with open(path, 'w') as f:
    json.dump(description, f, indent=4)
"
    fi

    echo "Starting $sim_name"
    (cd "$sim_dir" && "$radiseq_bin" parameters.txt > run.log 2>&1)  &
    pids+=($!)
    names+=("$sim_name")
done

if [ ${#pids[@]} -eq 0 ]; then
    echo "Nothing to run."
    exit 0
fi

exit_code=0
for i in "${!pids[@]}"; do
    if wait "${pids[$i]}"; then
        echo "Finished ${names[$i]}"
    else
        echo "FAILED ${names[$i]} (see $run_results_dir/${names[$i]}/run.log)"
        exit_code=1
    fi
done

exit $exit_code
