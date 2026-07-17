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
    expect_status 64 utility-selector env \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" UTILITY CnvMath
    expect_status 64 workflow-selector env \
        EPOLYSCAT_EXECUTABLE="$fake_executable" \
        bash "$wrapper" WORKFLOW ePolyScat_Run
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
