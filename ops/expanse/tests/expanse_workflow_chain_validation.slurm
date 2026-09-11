#!/bin/bash

#SBATCH --job-name=eps-workflow-chain
#SBATCH --account=nsf106
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=4
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

: "${CONTROLLER:?CONTROLLER must point to the formal ePolyScat controller}"
: "${OPENMOLCAS_INPUT:?OPENMOLCAS_INPUT must point to an OpenMolcas input}"
: "${EPOLYSCAT_INPUT:?EPOLYSCAT_INPUT must point to an ePolyScat input}"
: "${CNVLINFULL_CONTROL:?CNVLINFULL_CONTROL must point to CnvLinFull control data}"
[ -s "$OPENMOLCAS_INPUT" ]
[ -s "$EPOLYSCAT_INPUT" ]
[ -s "$CNVLINFULL_CONTROL" ]

validation_root="$SLURM_SUBMIT_DIR/results-workflow-chain"
stage1="$validation_root/01-data-generation"
stage2="$validation_root/02-epolyscat-run"
stage3="$validation_root/03-analysis"
rm -rf "$validation_root"
mkdir -p "$stage1" "$stage2" "$stage3"

cp "$OPENMOLCAS_INPUT" "$stage1/"
molcas_input=$(basename "$OPENMOLCAS_INPUT")
cd "$stage1"
"$CONTROLLER" WORKFLOW Data_Gen OpenMolcas "$molcas_input" \
    > controller.log 2>&1
grep -qiE '_RC_ALL_IS_WELL_|Happy landing' "$molcas_input.out"
molcas_molden=$(find . -maxdepth 1 -type f -name '*.scf.molden' -size +0c \
    -print -quit)
[ -n "$molcas_molden" ]
printf 'PASS Data_Gen OpenMolcas %s\n' "$(basename "$molcas_molden")" > PASS

cp "$EPOLYSCAT_INPUT" "$stage2/test12-with-dump.inp"
cp "$molcas_molden" "$stage2/test12.molcas.scf.molden"
cd "$stage2"
export pt=.
export po=.
export pd=.
"$CONTROLLER" WORKFLOW ePolyScat_Run test12-with-dump.inp \
    > controller.log 2>&1
grep -q 'Finalize' ePolyScat.out
grep -qE 'End (CrossSection|EDCS)' ePolyScat.out
[ -s test12.dump ]
grep -q 'DumpIdy Kinetic Energy' test12.dump
printf 'PASS ePolyScat_Run test12.dump\n' > PASS

cp "$CNVLINFULL_CONTROL" "$stage3/cnvlinfull.control"
cp test12.dump "$stage3/dumpidy.dat"
cd "$stage3"
"$CONTROLLER" WORKFLOW Analysis CnvLinFull \
    cnvlinfull.control dumpidy.dat > controller.log 2>&1
[ -s cnvlinfull.dat ]
grep -q 'Total cross section from coefficients' CnvLinFull.out
grep -q 'Total cross section from integration' CnvLinFull.out
printf 'PASS Analysis CnvLinFull cnvlinfull.dat\n' > PASS

{
    cat "$stage1/PASS"
    cat "$stage2/PASS"
    cat "$stage3/PASS"
    sha256sum \
        "$stage2/test12.molcas.scf.molden" \
        "$stage2/test12.dump" \
        "$stage3/dumpidy.dat" \
        "$stage3/cnvlinfull.dat"
} | tee "$validation_root/CHAIN-PASS"
