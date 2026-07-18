#!/usr/bin/env bash

set -euo pipefail

root=$(cd "$(dirname "$0")/.." && pwd)
wrapper="$root/epolyscat_wrapper.sh"

if [[ ! -f "$wrapper" ]]; then
    printf 'missing wrapper under test: %s\n' "$wrapper" >&2
    exit 1
fi

work=$(mktemp -d)
trap 'rm -rf "$work"' EXIT

fake_executable="$work/fake-ePolyScat.exe"
fake_launcher="$work/fake-ibrun"
utility_bin="$work/utility-bin"
mkdir -p "$utility_bin"

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
exec "$@"
EOF

chmod 700 "$fake_executable" "$fake_launcher"

for utility in CnvMath CnvMatLab CnvLinFull MoldenMerge NRFPAD Cube2igor; do
    cat > "$utility_bin/$utility.exe" <<'EOF'
#!/bin/sh

utility=$(basename "$0" .exe)
cat > "$UTILITY_LOG_DIR/$utility.stdin"
printf '%s\n' "$*" > "$UTILITY_LOG_DIR/$utility.args"

if [ "$utility" = "Cube2igor" ]; then
    [ -f fort.10 ] || {
        echo "fake Cube2igor: fort.10 was not staged" >&2
        exit 93
    }
    cmp -s fort.10 "$EXPECTED_CUBE_INPUT" || {
        echo "fake Cube2igor: fort.10 differs from the cube input" >&2
        exit 94
    }
    printf 'IGOR grid\n' > fort.11
    printf 'IGOR geometry\n' > fort.12
fi

printf '%s completed\n' "$utility"
EOF
    chmod 700 "$utility_bin/$utility.exe"
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
	EPOLYSCAT_EXECUTABLE="$fake_executable" \
    EPOLYSCAT_LAUNCHER="$fake_launcher" \
    LAUNCH_LOG="$work/portal-launch.log" \
        bash "$wrapper" MODULE ePolyScat
    grep -q '^Finalize$' ePolyScat.out
    [[ ! -e inp.dat ]]
    [[ ! -e injected ]]
    [[ "$(cat "$work/portal-launch.log")" == "$fake_executable" ]]
)

direct_dir="$work/direct"
mkdir -p "$direct_dir"
printf 'LMax 8\n\nGetBlms\n' > "$direct_dir/custom.inp"
cp "$direct_dir/custom.inp" "$work/direct-expected"
(
    cd "$direct_dir"
    EXPECTED_INPUT_FILE="$work/direct-expected" \
    EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" custom.inp
    grep -q '^End EDCS$' ePolyScat.out
    [[ ! -e inp.dat ]]
)

missing_dir="$work/missing"
mkdir -p "$missing_dir"
(
    cd "$missing_dir"
    expect_status 66 missing-input env \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" MODULE ePolyScat
)

unsupported_dir="$work/unsupported"
mkdir -p "$unsupported_dir"
printf 'LMax 4\n' > "$unsupported_dir/ep_input"
(
    cd "$unsupported_dir"
    expect_status 64 gaussian-selector env \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" MODULE Gaussian16
    expect_status 64 workflow-selector env \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" WORKFLOW ePolyScat_Run
)

utility_dir="$work/utilities"
utility_logs="$work/utility-logs"
mkdir -p "$utility_dir" "$utility_logs"
printf 'control line one\ncontrol line two\n' > "$utility_dir/control.in"
printf 'bend data\n' > "$utility_dir/bendorient.dat"
printf 'dump data\n' > "$utility_dir/dumpidy.dat"
printf 'first molden\n' > "$utility_dir/first.molden"
printf 'second molden\n' > "$utility_dir/second.molden"
printf 'orient data\n' > "$utility_dir/orientncro.dat"
printf 'cube data\n' > "$utility_dir/density.cube"

(
    cd "$utility_dir"
    for specification in \
        "CnvMath bendorient.dat" \
        "CnvMatLab bendorient.dat" \
        "CnvLinFull dumpidy.dat" \
        "NRFPAD orientncro.dat"; do
        set -- $specification
        utility=$1
        data_file=$2
        UTILITY_LOG_DIR="$utility_logs" \
        EPOLYSCAT_UTILITY_BIN_DIR="$utility_bin" \
            bash "$wrapper" UTILITY "$utility" control.in "$data_file"
        cmp -s control.in "$utility_logs/$utility.stdin"
        grep -q "^$utility completed$" "$utility.out"
    done

    UTILITY_LOG_DIR="$utility_logs" \
    EPOLYSCAT_UTILITY_BIN_DIR="$utility_bin" \
        bash "$wrapper" UTILITY MoldenMerge control.in first.molden second.molden
    cmp -s control.in "$utility_logs/MoldenMerge.stdin"

    EXPECTED_CUBE_INPUT="$utility_dir/density.cube" \
    UTILITY_LOG_DIR="$utility_logs" \
    EPOLYSCAT_UTILITY_BIN_DIR="$utility_bin" \
        bash "$wrapper" UTILITY Cube2igor control.in density.cube
    cmp -s control.in "$utility_logs/Cube2igor.stdin"
    grep -q '^IGOR grid$' fort.11
    grep -q '^IGOR geometry$' fort.12
    [[ ! -e fort.10 ]]

    UTILITY_LOG_DIR="$utility_logs" \
    EPOLYSCAT_UTILITY_BIN_DIR="$utility_bin" \
        bash "$wrapper" WORKFLOW Analysis CnvLinFull control.in dumpidy.dat

    expect_status 64 utility-missing-control env \
        EPOLYSCAT_UTILITY_BIN_DIR="$utility_bin" \
        bash "$wrapper" UTILITY CnvMath
    expect_status 64 molden-one-input env \
        EPOLYSCAT_UTILITY_BIN_DIR="$utility_bin" \
        bash "$wrapper" UTILITY MoldenMerge control.in first.molden
)

occupied_dir="$work/occupied"
mkdir -p "$occupied_dir"
printf 'LMax 4\n' > "$occupied_dir/ep_input"
printf 'do not overwrite\n' > "$occupied_dir/inp.dat"
(
    cd "$occupied_dir"
    expect_status 73 occupied-input env \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" MODULE ePolyScat
    grep -q '^do not overwrite$' inp.dat
)

failure_dir="$work/failures"
mkdir -p "$failure_dir"
printf 'LMax 4\n' > "$failure_dir/ep_input"
(
    cd "$failure_dir"
    expect_status 42 runner-failure env \
        FAKE_MODE=fail EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" MODULE ePolyScat
    [[ ! -e inp.dat ]]

    expect_status 70 abnormal-output env \
        FAKE_MODE=abnormal EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" MODULE ePolyScat
    [[ ! -e inp.dat ]]

    expect_status 70 incomplete-output env \
        FAKE_MODE=incomplete EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" MODULE ePolyScat
    [[ ! -e inp.dat ]]
)

printf 'Frontera ePolyScat wrapper regression tests passed\n'
