#!/bin/bash

#SBATCH --job-name=eps-cnvlinfull
#SBATCH --account=nsf106
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:10:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

: "${CONTROLLER:?CONTROLLER must point to the staged ePolyScat controller}"
: "${CONTROL_FILE:?CONTROL_FILE must point to the CnvLinFull control file}"
: "${DUMP_FILE:?DUMP_FILE must point to real DumpIdy output}"

validation_root="$SLURM_SUBMIT_DIR/results-cnvlinfull"
rm -rf "$validation_root"
mkdir -p "$validation_root"

run_case() {
    route=$1
    case_dir="$validation_root/$route"
    mkdir -p "$case_dir"
    cp "$CONTROL_FILE" "$case_dir/cnvlinfull.control"
    cp "$DUMP_FILE" "$case_dir/dumpidy.dat"
    cd "$case_dir"

    if [ "$route" = "UTILITY" ]; then
        "$CONTROLLER" UTILITY CnvLinFull cnvlinfull.control dumpidy.dat \
            > controller.log 2>&1
    else
        "$CONTROLLER" WORKFLOW Analysis CnvLinFull \
            cnvlinfull.control dumpidy.dat > controller.log 2>&1
    fi

    [ -s cnvlinfull.dat ]
    grep -q 'Total cross section from coefficients' CnvLinFull.out
    grep -q 'Total cross section from integration' CnvLinFull.out
    printf 'PASS %s CnvLinFull\n' "$route" | tee PASS
}

run_case UTILITY
run_case WORKFLOW
