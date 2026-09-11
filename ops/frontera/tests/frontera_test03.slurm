#!/bin/bash
#SBATCH -J eps-wrapper-test03
#SBATCH -p development
#SBATCH -N 1
#SBATCH -n 4
#SBATCH -A PHY20028
#SBATCH -t 00:30:00
#SBATCH -o ePolyScat.stdout
#SBATCH -e ePolyScat.stderr

set -u

cd "$SLURM_SUBMIT_DIR"

export po="$SLURM_SUBMIT_DIR"
export pt="$SLURM_SUBMIT_DIR"
export pd="$SLURM_SUBMIT_DIR"

exec /home1/06170/ampuser/apps/ePolyScat/ePolyScat MODULE ePolyScat
