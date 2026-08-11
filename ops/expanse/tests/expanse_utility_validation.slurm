#!/bin/bash

#SBATCH --job-name=eps-utility-matrix
#SBATCH --account=nsf106
#SBATCH --partition=shared
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --time=00:15:00
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

: "${CONTROLLER:?CONTROLLER must point to the staged ePolyScat controller}"
: "${FIXTURES:?FIXTURES must point to the validation fixture directory}"
: "${MOLDEN_FILE:?MOLDEN_FILE must point to a complete Molden file}"

validation_root="$SLURM_SUBMIT_DIR/results"
mkdir -p "$validation_root"

require_nonempty() {
    [ -s "$1" ] || {
        echo "Expected non-empty output: $1" >&2
        return 1
    }
}

require_matching_output() {
    pattern=$1
    find . -maxdepth 1 -type f -name "$pattern" -size +0c | grep -q . || {
        echo "Expected a non-empty output matching: $pattern" >&2
        return 1
    }
}

prepare_case() {
    case_dir=$1
    rm -rf "$case_dir"
    mkdir -p "$case_dir"
}

invoke_controller() {
    route=$1
    utility=$2
    shift 2
    if [ "$route" = "UTILITY" ]; then
        "$CONTROLLER" UTILITY "$utility" "$@"
    else
        "$CONTROLLER" WORKFLOW Analysis "$utility" "$@"
    fi
}

run_utility_case() {
    route=$1
    utility=$2
    case_dir="$validation_root/${route}-${utility}"
    prepare_case "$case_dir"
    cd "$case_dir"

    case "$utility" in
    CnvMath)
        cp "$FIXTURES/cnvmath.control" .
        cp "$FIXTURES/orient-coefficients.dat" .
        invoke_controller "$route" CnvMath cnvmath.control orient-coefficients.dat \
            > controller.log 2>&1
        require_matching_output '*.m'
        require_matching_output '*.tex'
        require_matching_output '*FLN.dat'
        grep -q 'Expansion in terms of Legendre Polynomials' CnvMath.out
        ;;
    CnvMatLab)
        cp "$FIXTURES/cnvmatlab.control" .
        cp "$FIXTURES/orient-coefficients.dat" .
        invoke_controller "$route" CnvMatLab cnvmatlab.control orient-coefficients.dat \
            > controller.log 2>&1
        require_matching_output '*L*.dat'
        require_matching_output '*FLN.dat'
        grep -q 'Expansion in terms of Legendre Polynomials' CnvMatLab.out
        ;;
    MoldenMerge)
        cp "$FIXTURES/moldenmerge.control" .
        cp "$MOLDEN_FILE" molden-a.molden
        cp "$MOLDEN_FILE" molden-b.molden
        invoke_controller "$route" MoldenMerge moldenmerge.control \
            molden-a.molden molden-b.molden > controller.log 2>&1
        require_nonempty merged.molden
        grep -q '^\[MO\]$' merged.molden
        grep -q 'Combine orbitals from different molden files' MoldenMerge.out
        ;;
    NRFPAD)
        cp "$FIXTURES/nrfpad.control" .
        cp "$FIXTURES/nrfpad-fln.dat" .
        invoke_controller "$route" NRFPAD nrfpad.control nrfpad-fln.dat \
            > controller.log 2>&1
        require_nonempty fort.30
        grep -q '^FLNNu' fort.30
        grep -q 'Energy Use' NRFPAD.out
        ;;
    Cube2igor)
        cp "$FIXTURES/cube2igor.control" .
        cp "$FIXTURES/validation.cube" .
        invoke_controller "$route" Cube2igor cube2igor.control validation.cube \
            > controller.log 2>&1
        require_nonempty fort.11
        require_nonempty fort.12
        grep -q '^Macro ValidationCube' fort.11
        grep -q '^IGOR$' fort.12
        grep -q 'Finshed reading orbital' Cube2igor.out
        ;;
    *)
        echo "Unsupported validation utility: $utility" >&2
        return 64
        ;;
    esac

    printf 'PASS %s %s\n' "$route" "$utility" | tee PASS
}

for utility in CnvMath CnvMatLab MoldenMerge NRFPAD Cube2igor; do
    run_utility_case UTILITY "$utility"
done

for utility in CnvMath CnvMatLab MoldenMerge NRFPAD Cube2igor; do
    run_utility_case WORKFLOW "$utility"
done

find "$validation_root" -name PASS -type f -print | sort
