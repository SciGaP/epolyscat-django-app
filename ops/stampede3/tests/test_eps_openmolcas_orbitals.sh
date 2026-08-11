#!/bin/bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
default_wrapper=$(cd -- "$script_dir/.." && pwd)/eps_openmolcas.sh
wrapper=${OPENMOLCAS_WRAPPER:-$default_wrapper}

fail() {
    echo "FAIL: $1" >&2
    exit 1
}

[[ -x "$wrapper" ]] || fail "OpenMolcas wrapper is missing or not executable: $wrapper"

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT
mkdir -p "$tmp_dir/bin" "$tmp_dir/run-success" "$tmp_dir/run-render-failure"

cat > "$tmp_dir/minimal.molden" <<'MOLDEN'
[Molden Format]
[Atoms] Angs
H 1 1 0.0 0.0 0.0
[GTO]
1 0
s 1 1.00
1.0 1.0
[MO]
Sym= 1
Ene= -0.5
Spin= Alpha
Occup= 2.0
1 1.0
MOLDEN

cat > "$tmp_dir/bin/pymolcas" <<'FAKE_MOLCAS'
#!/bin/bash
set -euo pipefail
input=$1
project=$(basename "$input")
project=${project%.*}
cp "$MOLDEN_FIXTURE" "$project.scf.molden"
printf '%s\n' 'Happy landing' > molcas.status
echo "fake OpenMolcas completed"
FAKE_MOLCAS

cat > "$tmp_dir/bin/renderer" <<'FAKE_RENDERER'
#!/bin/bash
set -euo pipefail
printf '%s\n' "$1" "$2" > "$RENDER_CALL_FILE"
exit "${RENDER_EXIT_CODE:-0}"
FAKE_RENDERER
chmod +x "$tmp_dir/bin/pymolcas" "$tmp_dir/bin/renderer"

run_wrapper() {
    local run_dir=$1
    local render_status=$2
    printf '%s\n' '&GATEWAY' > "$run_dir/molcas.inp"
    (
        cd "$run_dir"
        PYMOLCAS="$tmp_dir/bin/pymolcas" \
        MOLDEN_FIXTURE="$tmp_dir/minimal.molden" \
        ORBITAL_RENDERER="$tmp_dir/bin/renderer" \
        RENDER_CALL_FILE="$run_dir/render-call" \
        RENDER_EXIT_CODE="$render_status" \
        MOLCAS_WORKDIR="$run_dir/tmp" \
            "$wrapper" molcas.inp > wrapper.stdout 2> wrapper.stderr
    )
}

run_wrapper "$tmp_dir/run-success" 0
[[ -s "$tmp_dir/run-success/molcas.inp.out" ]] || fail "OpenMolcas output log was not created"
[[ -s "$tmp_dir/run-success/render-call" ]] || fail "orbital renderer was not called"
grep -qx 'molcas.scf.molden' "$tmp_dir/run-success/render-call" ||
    fail "renderer did not receive the SCF Molden file"
grep -qx "$tmp_dir/run-success" "$tmp_dir/run-success/render-call" ||
    fail "renderer did not receive the run output directory"

run_wrapper "$tmp_dir/run-render-failure" 70
[[ -s "$tmp_dir/run-render-failure/molcas.inp.out" ]] ||
    fail "scientific output was lost when optional rendering failed"
grep -q 'optional orbital image generation failed; calculation outputs are preserved' \
    "$tmp_dir/run-render-failure/wrapper.stderr" ||
    fail "non-fatal renderer failure warning was not emitted"

echo "PASS: eps_openmolcas orbital integration"
