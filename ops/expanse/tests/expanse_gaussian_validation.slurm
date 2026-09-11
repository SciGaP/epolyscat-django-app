#!/bin/bash

#SBATCH --job-name=eps-gaussian-check
#SBATCH --account=nsf106
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mem=2G
#SBATCH --time=00:05:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -u

: "${CONTROLLER:?CONTROLLER must point to the staged ePolyScat controller}"
: "${GAUSSIAN_INPUT:?GAUSSIAN_INPUT must point to a Gaussian input file}"
: "${VALIDATION_ROUTE:?VALIDATION_ROUTE must be MODULE or WORKFLOW}"

work_dir="$SLURM_SUBMIT_DIR/results-gaussian-$VALIDATION_ROUTE"
rm -rf "$work_dir"
mkdir -p "$work_dir"
cp "$GAUSSIAN_INPUT" "$work_dir/gaussian.in"
cd "$work_dir"

module --force purge > module-preflight.log 2>&1
module load cpu/0.15.4 gaussian/16.C.01 >> module-preflight.log 2>&1
module_status=$?
printf 'Gaussian module load status: %s\n' "$module_status" \
    >> module-preflight.log

export AIRAVATA_USERNAME=zzhang90
export AIRAVATA_INPUTS='[gaussian.in]'
export AIRAVATA_ExptDataDir=

if [ "$VALIDATION_ROUTE" = "MODULE" ]; then
    "$CONTROLLER" MODULE Gaussian16 gaussian.in > controller.log 2>&1
    controller_status=$?
else
    "$CONTROLLER" WORKFLOW Data_Gen Gaussian16 gaussian.in \
        > controller.log 2>&1
    controller_status=$?
fi

printf 'Controller status: %s\n' "$controller_status" | tee RESULT
if [ "$module_status" -eq 0 ] && [ "$controller_status" -eq 0 ]; then
    grep -q 'Normal termination of Gaussian' gaussian.log
    printf 'PASS %s Gaussian16\n' "$VALIDATION_ROUTE" | tee PASS
    exit 0
fi

printf 'BLOCKED %s Gaussian16 module=%s controller=%s\n' \
    "$VALIDATION_ROUTE" "$module_status" "$controller_status" | tee BLOCKED
exit "${controller_status:-69}"
