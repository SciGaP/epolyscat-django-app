#!/bin/bash

#SBATCH --job-name=eps-openmolcas
#SBATCH --account=nsf106
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --time=00:20:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

: "${CONTROLLER:?CONTROLLER must point to the staged ePolyScat controller}"
: "${MOLCAS_INPUT:?MOLCAS_INPUT must point to an OpenMolcas input file}"
: "${VALIDATION_ROUTE:?VALIDATION_ROUTE must be MODULE or WORKFLOW}"

case "$VALIDATION_ROUTE" in
MODULE|WORKFLOW)
    ;;
*)
    echo "Unsupported validation route: $VALIDATION_ROUTE" >&2
    exit 64
    ;;
esac

work_dir="$SLURM_SUBMIT_DIR/results-openmolcas-$VALIDATION_ROUTE"
rm -rf "$work_dir"
mkdir -p "$work_dir"
cp "$MOLCAS_INPUT" "$work_dir/"
cd "$work_dir"

input_name=$(basename "$MOLCAS_INPUT")
if [ "$VALIDATION_ROUTE" = "MODULE" ]; then
    "$CONTROLLER" MODULE OpenMolcas "$input_name" > controller.log 2>&1
else
    "$CONTROLLER" WORKFLOW Data_Gen OpenMolcas "$input_name" \
        > controller.log 2>&1
fi

grep -qiE '_RC_ALL_IS_WELL_|Happy landing' "$input_name.out"
find . -maxdepth 1 -type f -name '*.scf.molden' -size +0c | grep -q .
find . -maxdepth 1 -type f -name '*.ScfOrb' -size +0c | grep -q .

printf 'PASS %s OpenMolcas\n' "$VALIDATION_ROUTE" | tee PASS
