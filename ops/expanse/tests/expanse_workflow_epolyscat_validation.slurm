#!/bin/bash

#SBATCH --job-name=eps-workflow-test12
#SBATCH --account=nsf106
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:20:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

: "${CONTROLLER:?CONTROLLER must point to the staged ePolyScat controller}"
: "${EPOLYSCAT_INPUT:?EPOLYSCAT_INPUT must point to the test12 input file}"
: "${EPOLYSCAT_MOLDEN:?EPOLYSCAT_MOLDEN must point to the test12 Molden file}"

work_dir="$SLURM_SUBMIT_DIR/results-workflow-epolyscat"
rm -rf "$work_dir"
mkdir -p "$work_dir"
cp "$EPOLYSCAT_INPUT" "$work_dir/test12-with-dump.inp"
cp "$EPOLYSCAT_MOLDEN" "$work_dir/test12.molcas.scf.molden"
cd "$work_dir"

printf '%s\n' \
    "FileName 'DumpOut' 'test12.dump' 'REWIND'" \
    "DumpIdy 'test12SU.idy' 10.0" \
    >> test12-with-dump.inp

export pt=.
export po=.
export pd=.

"$CONTROLLER" WORKFLOW ePolyScat_Run test12-with-dump.inp \
    > controller.log 2>&1

grep -q 'Finalize' ePolyScat.out
[ -s test12SU.idy ]
[ -s test12PU.idy ]
[ -s test12.dump ]
grep -q 'DumpIdy Kinetic Energy' test12.dump

printf 'PASS WORKFLOW ePolyScat_Run\n' | tee PASS
