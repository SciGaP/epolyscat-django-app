#!/bin/bash

set -uo pipefail

fail() {
    echo "OpenMolcas wrapper: $1" >&2
    exit "${2:-1}"
}

input=${1-}
[ -n "$input" ] || fail "usage: $0 molcas_input" 64
[ -r "$input" ] || fail "input file is missing or unreadable: $input" 66

nprocs=${SLURM_NPROCS:-${SLURM_NTASKS:-1}}
export MOLCAS_NPROCS=$nprocs
export OMP_NUM_THREADS=${SLURM_NTASKS:-$nprocs}
export OMP_STACKSIZE=${OMP_STACKSIZE:-140736864318339}
export OMP_PROC_BIND=${OMP_PROC_BIND:-spread}
ulimit -s unlimited 2>/dev/null || true

job_id=${SLURM_JOB_ID:-${SLURM_JOBID:-$$}}
scratch_candidate="/scratch/${USER:-ampuser}/job_$job_id"
if [[ -d "$scratch_candidate" ]]; then
    export MOLCAS_WORKDIR=$scratch_candidate
else
    export MOLCAS_WORKDIR=${MOLCAS_WORKDIR:-$PWD/tmp}
    mkdir -p "$MOLCAS_WORKDIR" ||
        fail "cannot create work directory: $MOLCAS_WORKDIR" 73
fi

if [[ -n "${SLURM_MEM_PER_NODE:-}" && "${SLURM_MEM_PER_NODE}" != "0" ]]; then
    export MOLCAS_MEM=$SLURM_MEM_PER_NODE
elif [[ -z "${MOLCAS_MEM:-}" || "${MOLCAS_MEM}" == "0" ]]; then
    unset MOLCAS_MEM
fi

if [[ -z "${PYMOLCAS:-}" ]]; then
    pyxchem_env=${PYXCHEM_ENV:-/home/ampuser/apps/xchem/pyxchem_new.env}
    if [[ -r "$pyxchem_env" ]]; then
        # shellcheck disable=SC1090
        source "$pyxchem_env"
    fi
fi

molcas_root=${MOLCAS_ROOT:-/home/ampuser/apps/xchem/OpenMolcas}
pymolcas=${PYMOLCAS:-$molcas_root/sbin/pymolcas}
[[ -x "$pymolcas" ]] ||
    fail "pymolcas is missing or not executable: $pymolcas" 69

output="${input}.out"
"$pymolcas" "$input" | tee "$output"
molcas_status=${PIPESTATUS[0]}
if [[ "$molcas_status" -ne 0 ]]; then
    fail "pymolcas exited with status $molcas_status; inspect $output for details" \
        "$molcas_status"
fi

if [[ -f molcas.status ]]; then
    grep -qi "Happy landing" molcas.status ||
        fail "molcas.status does not report Happy landing" 70
elif ! grep -qiE "_RC_ALL_IS_WELL_|Happy landing" "$output"; then
    fail "OpenMolcas output does not contain a successful completion marker" 70
fi

project=$(basename "$input")
project=${project%.*}
molden_file=""
if [[ -s "$project.rasscf.molden" ]]; then
    molden_file="$project.rasscf.molden"
elif [[ -s "$project.scf.molden" ]]; then
    molden_file="$project.scf.molden"
else
    fail "OpenMolcas did not produce a non-empty Molden orbital file" 70
fi

render_orbitals=${3:-yes}
if [[ "$render_orbitals" != "no" ]]; then
    wrapper_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
    orbital_renderer=${ORBITAL_RENDERER:-$wrapper_dir/render_molden_orbitals.sh}
    if [[ ! -x "$orbital_renderer" ]]; then
        echo "OpenMolcas wrapper: orbital image generation skipped because the optional renderer is unavailable" >&2
    elif ! "$orbital_renderer" "$molden_file" "$PWD"; then
        echo "OpenMolcas wrapper: optional orbital image generation failed; calculation outputs are preserved" >&2
    fi
fi

exit 0
