#!/bin/bash

#SBATCH --job-name=eps-controller-test03
#SBATCH --account=nsf106
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -u

: "${CONTROLLER:?CONTROLLER must point to the staged ePolyScat controller}"

cd "$SLURM_SUBMIT_DIR"
export pt=.
export po=.
export pd=.

"$CONTROLLER" MODULE ePolyScat
