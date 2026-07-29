#!/usr/bin/env bash

set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
controller="$root/epolyscat_controller.sh"

if [[ ! -f "$controller" ]]; then
    printf 'missing controller under test: %s\n' "$controller" >&2
    exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

fake_executable="$work/fake-ePolyScat.exe"
fake_launcher="$work/fake-mpirun"
fake_bin="$work/fake-bin"
fake_runtime_dir="$fake_bin/expanseintelmkl"
delegate_log="$work/delegate.log"
utility_log_dir="$work/utility-logs"
mkdir -p "$fake_runtime_dir" "$utility_log_dir"

cat > "$fake_executable" <<'EOF'
#!/bin/sh

if [ -n "${EXPECTED_INPUT_FILE-}" ] && ! cmp -s inp.dat "$EXPECTED_INPUT_FILE"; then
    echo "fake ePolyScat: inp.dat differs from the staged input" >&2
    exit 91
fi

case "${FAKE_MODE-ok}" in
ok)
    printf 'End EDCS\nFinalize\n'
    ;;
fail)
    printf 'runner failed\n'
    exit 42
    ;;
abnormal)
    printf 'Abnormal Ending\n'
    ;;
incomplete)
    printf 'calculation stopped before finalization\n'
    ;;
*)
    echo "unknown fake mode: ${FAKE_MODE}" >&2
    exit 92
    ;;
esac
EOF

cat > "$fake_launcher" <<'EOF'
#!/bin/sh

printf '%s\n' "$*" > "$LAUNCH_LOG"
while [ "${1-}" = "-bootstrap" ] || [ "${1-}" = "-np" ]; do
    shift 2
done
exec "$@"
EOF

cat > "$fake_bin/eps_g16.sh" <<'EOF'
#!/bin/sh
printf 'Gaussian16:%s\n' "$*" >> "$DELEGATE_LOG"
case "${FAKE_DELEGATE_MODE-ok}" in
ok)
    output_file=${1%.*}.log
    printf 'Normal termination of Gaussian 16\n' > "$output_file"
    ;;
false-success)
    exit 0
    ;;
*)
    exit 41
    ;;
esac
EOF

cat > "$fake_bin/eps_openmolcas.sh" <<'EOF'
#!/bin/sh
printf 'OpenMolcas:%s\n' "$*" >> "$DELEGATE_LOG"
[ "${FAKE_DELEGATE_MODE-ok}" = "ok" ] || exit 43
EOF

chmod 700 \
    "$fake_executable" \
    "$fake_launcher" \
    "$fake_bin/eps_g16.sh" \
    "$fake_bin/eps_openmolcas.sh"

for utility in CnvMath CnvMatLab CnvLinFull MoldenMerge NRFPAD Cube2igor; do
    cat > "$fake_runtime_dir/$utility.exe" <<'EOF'
#!/bin/sh

utility=$(basename "$0" .exe)
cat > "$UTILITY_LOG_DIR/$utility.stdin"

if [ "$utility" = "Cube2igor" ]; then
    [ -f fort.10 ] || exit 93
    cmp -s fort.10 "$EXPECTED_CUBE_INPUT" || exit 94
    printf 'IGOR grid\n' > fort.11
    printf 'IGOR geometry\n' > fort.12
fi

if [ "$utility" = "NRFPAD" ]; then
    [ -f fort.21 ] || exit 95
    cmp -s fort.21 "$EXPECTED_NRFPAD_FLN" || exit 96
    if [ -n "${EXPECTED_NRFPAD_HNU-}" ]; then
        [ -f fort.22 ] || exit 97
        cmp -s fort.22 "$EXPECTED_NRFPAD_HNU" || exit 98
    fi
    printf 'RFPAD output\n' > fort.20
    printf 'Matlab coefficients\n' > fort.23
fi

if [ "${FAKE_UTILITY_MODE-ok}" = "abstop" ]; then
    printf 'ABSTOP: invalid scientific input\n'
    exit 0
fi

printf '%s completed\n' "$utility"
EOF
    chmod 700 "$fake_runtime_dir/$utility.exe"
done

expect_status() {
    local expected=$1
    local label=$2
    shift 2

    set +e
    "$@" > "$work/$label.stdout" 2> "$work/$label.stderr"
    local actual=$?
    set -e

    if [[ "$actual" -ne "$expected" ]]; then
        printf '%s: expected status %s, got %s\n' "$label" "$expected" "$actual" >&2
        cat "$work/$label.stdout" >&2
        cat "$work/$label.stderr" >&2
        exit 1
    fi
}

portal_dir="$work/portal"
mkdir -p "$portal_dir"
cat > "$portal_dir/ep_input" <<'EOF'
LMax 15

ScatEng 3.0 4.0
Convert '$pt/test03.molden2012' 'molden'
DumpIdy '$po/scattering.idy'
SaveData '$pd/scattering.dat'
# Literal shell text must not execute: $(touch injected)


EOF
cat > "$work/portal-expected" <<'EOF'
LMax 15

ScatEng 3.0 4.0
Convert './test03.molden2012' 'molden'
DumpIdy './scattering.idy'
SaveData './scattering.dat'
# Literal shell text must not execute: $(touch injected)


EOF

(
    cd "$portal_dir"
    EXPECTED_INPUT_FILE="$work/portal-expected" \
    pt=. po=. pd=. \
    SLURM_JOB_ID=123 SLURM_NTASKS=4 \
    EPOLYSCAT_SKIP_MODULE_SETUP=1 \
    EPOLYSCAT_EXECUTABLE="$fake_executable" \
    EPOLYSCAT_LAUNCHER="$fake_launcher" \
    EPOLYSCAT_BIN_DIR="$fake_bin" \
    LAUNCH_LOG="$work/portal-launch.log" \
        sh "$controller" MODULE ePolyScat
    grep -q '^Finalize$' ePolyScat.out
    [[ ! -e inp.dat ]]
    [[ ! -e injected ]]
    [[ "$(cat "$work/portal-launch.log")" == "-bootstrap slurm -np 4 $fake_executable" ]]
)

workflow_dir="$work/workflow"
mkdir -p "$workflow_dir"
printf 'LMax 8\n\nGetBlms\n' > "$workflow_dir/ep_input"
cp "$workflow_dir/ep_input" "$work/workflow-expected"
(
    cd "$workflow_dir"
    EXPECTED_INPUT_FILE="$work/workflow-expected" \
    EPOLYSCAT_SKIP_MODULE_SETUP=1 \
    EPOLYSCAT_EXECUTABLE="$fake_executable" \
    EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" WORKFLOW ePolyScat_Run
    grep -q '^End EDCS$' ePolyScat.out
    [[ ! -e inp.dat ]]
)

delegate_dir="$work/delegates"
mkdir -p "$delegate_dir"
printf 'gaussian input\n' > "$delegate_dir/gaussian.in"
printf 'molcas input\n' > "$delegate_dir/molcas.inp"
(
    cd "$delegate_dir"
    DELEGATE_LOG="$delegate_log" EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" MODULE Gaussian16 gaussian.in
    DELEGATE_LOG="$delegate_log" EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" WORKFLOW Data_Gen Gaussian16 gaussian.in
    DELEGATE_LOG="$delegate_log" EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" MODULE OpenMolcas molcas.inp
    DELEGATE_LOG="$delegate_log" EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" WORKFLOW Data_Gen OpenMolcas molcas.inp
)
[[ "$(grep -c '^Gaussian16:gaussian.in$' "$delegate_log")" -eq 2 ]]
[[ "$(grep -c '^OpenMolcas:molcas.inp -i=yes$' "$delegate_log")" -eq 2 ]]

(
    cd "$delegate_dir"
    expect_status 41 gaussian-failure env \
        FAKE_DELEGATE_MODE=fail \
        DELEGATE_LOG="$delegate_log" \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" MODULE Gaussian16 gaussian.in
    rm -f gaussian.log
    expect_status 70 gaussian-false-success env \
        FAKE_DELEGATE_MODE=false-success \
        DELEGATE_LOG="$delegate_log" \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" MODULE Gaussian16 gaussian.in
)

utility_dir="$work/utilities"
mkdir -p "$utility_dir"
printf 'control line one\ncontrol line two\n' > "$utility_dir/control.in"
printf 'bend data\n' > "$utility_dir/bendorient.dat"
printf 'dump data\n' > "$utility_dir/dumpidy.dat"
printf 'first molden\n' > "$utility_dir/first.molden"
printf 'second molden\n' > "$utility_dir/second.molden"
printf 'orient data\n' > "$utility_dir/orientncro.dat"
printf 'dipole data\n' > "$utility_dir/hnu.dat"
printf 'cube data\n' > "$utility_dir/density.cube"

(
    cd "$utility_dir"
    for specification in \
        "CnvMath bendorient.dat" \
        "CnvMatLab bendorient.dat" \
        "CnvLinFull dumpidy.dat"; do
        set -- $specification
        utility=$1
        data_file=$2
        UTILITY_LOG_DIR="$utility_log_dir" \
        EPOLYSCAT_SKIP_MODULE_SETUP=1 \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
            sh "$controller" UTILITY "$utility" control.in "$data_file"
        cmp -s control.in "$utility_log_dir/$utility.stdin"
        grep -q "^$utility completed$" "$utility.out"
    done

    UTILITY_LOG_DIR="$utility_log_dir" \
    EXPECTED_NRFPAD_FLN="$utility_dir/orientncro.dat" \
    EXPECTED_NRFPAD_HNU="$utility_dir/hnu.dat" \
    EPOLYSCAT_SKIP_MODULE_SETUP=1 \
    EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" UTILITY NRFPAD control.in orientncro.dat hnu.dat
    grep -q '^RFPAD output$' fort.20
    [[ ! -e fort.21 ]]
    [[ ! -e fort.22 ]]

    UTILITY_LOG_DIR="$utility_log_dir" \
    EPOLYSCAT_SKIP_MODULE_SETUP=1 \
    EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" UTILITY MoldenMerge control.in first.molden second.molden

    UTILITY_LOG_DIR="$utility_log_dir" \
    EXPECTED_CUBE_INPUT="$utility_dir/density.cube" \
    EPOLYSCAT_SKIP_MODULE_SETUP=1 \
    EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" UTILITY Cube2igor control.in density.cube
    grep -q '^IGOR grid$' fort.11
    [[ ! -e fort.10 ]]

    UTILITY_LOG_DIR="$utility_log_dir" \
    EPOLYSCAT_SKIP_MODULE_SETUP=1 \
    EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" WORKFLOW Analysis CnvLinFull control.in dumpidy.dat

    expect_status 64 utility-missing-control env \
        EPOLYSCAT_SKIP_MODULE_SETUP=1 \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" UTILITY CnvMath
    expect_status 64 molden-one-input env \
        EPOLYSCAT_SKIP_MODULE_SETUP=1 \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" UTILITY MoldenMerge control.in first.molden
    expect_status 70 utility-abstop env \
        FAKE_UTILITY_MODE=abstop \
        UTILITY_LOG_DIR="$utility_log_dir" \
        EPOLYSCAT_SKIP_MODULE_SETUP=1 \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" UTILITY CnvLinFull control.in dumpidy.dat
)

failure_dir="$work/failures"
mkdir -p "$failure_dir"
printf 'LMax 4\n' > "$failure_dir/ep_input"
(
    cd "$failure_dir"
    expect_status 42 runner-failure env \
        FAKE_MODE=fail \
        EPOLYSCAT_SKIP_MODULE_SETUP=1 \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" MODULE ePolyScat
    [[ ! -e inp.dat ]]

    expect_status 70 abnormal-output env \
        FAKE_MODE=abnormal \
        EPOLYSCAT_SKIP_MODULE_SETUP=1 \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" MODULE ePolyScat

    expect_status 70 incomplete-output env \
        FAKE_MODE=incomplete \
        EPOLYSCAT_SKIP_MODULE_SETUP=1 \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        EPOLYSCAT_BIN_DIR="$fake_bin" \
        sh "$controller" MODULE ePolyScat
)

printf 'Expanse ePolyScat controller regression tests passed\n'
