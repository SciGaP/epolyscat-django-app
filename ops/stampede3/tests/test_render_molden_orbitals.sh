#!/bin/bash

set -euo pipefail

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
ops_dir=$(cd -- "$script_dir/.." && pwd)
renderer="$ops_dir/render_molden_orbitals.sh"
template="$ops_dir/print_frontier_orbitals.template"

fail() {
    echo "FAIL: $1" >&2
    exit 1
}

[[ -x "$renderer" ]] || fail "renderer is missing or not executable: $renderer"
[[ -r "$template" ]] || fail "template is missing or unreadable: $template"

tmp_dir=$(mktemp -d)
trap 'rm -rf "$tmp_dir"' EXIT

mkdir -p "$tmp_dir/fake-bin" "$tmp_dir/supported-output" "$tmp_dir/unsupported-output"
: > "$tmp_dir/fake-jmol.jar"
: > "$tmp_dir/java-calls"

cat > "$tmp_dir/fake-bin/java" <<'FAKE_JAVA'
#!/bin/bash
set -euo pipefail
script=${!#}
printf '%s\n' "$script" >> "$JMOL_TEST_CALLS"
grep -q 'FILTER "NOSORT"' "$script"
awk -F'"' '/write PNG/ { print $2 }' "$script" | while IFS= read -r output; do
    printf '\211PNG\r\n\032\n' > "$output"
done
FAKE_JAVA
chmod +x "$tmp_dir/fake-bin/java"

cat > "$tmp_dir/supported.molden" <<'MOLDEN'
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

cat > "$tmp_dir/unsupported.molden" <<'MOLDEN'
[Molden Format]
[9G]
[Atoms] Angs
H 1 1 0.0 0.0 0.0
MOLDEN

PATH="$tmp_dir/fake-bin:$PATH" \
JMOL_JAR="$tmp_dir/fake-jmol.jar" \
JMOL_TEMPLATE="$template" \
JMOL_TEST_CALLS="$tmp_dir/java-calls" \
    "$renderer" "$tmp_dir/supported.molden" "$tmp_dir/supported-output"

expected_images=(
    orbital_homo.png
    orbital_lumo.png
    orbital_homo_x.png
    orbital_lumo_x.png
    orbital_homo_y.png
    orbital_lumo_y.png
    orbital_homo_z.png
    orbital_lumo_z.png
)
for image in "${expected_images[@]}"; do
    [[ -s "$tmp_dir/supported-output/$image" ]] || fail "missing supported image: $image"
done

[[ $(wc -l < "$tmp_dir/java-calls") -eq 4 ]] || fail "expected four Jmol invocations"

before_calls=$(wc -l < "$tmp_dir/java-calls")
PATH="$tmp_dir/fake-bin:$PATH" \
JMOL_JAR="$tmp_dir/fake-jmol.jar" \
JMOL_TEMPLATE="$template" \
JMOL_TEST_CALLS="$tmp_dir/java-calls" \
    "$renderer" "$tmp_dir/unsupported.molden" "$tmp_dir/unsupported-output" \
    > "$tmp_dir/unsupported.stdout" 2> "$tmp_dir/unsupported.stderr"
after_calls=$(wc -l < "$tmp_dir/java-calls")

[[ "$before_calls" -eq "$after_calls" ]] || fail "Jmol ran for an unsupported basis"
grep -q 'unsupported angular momentum marker \[9G\]' "$tmp_dir/unsupported.stderr" ||
    fail "unsupported basis warning was not emitted"
if find "$tmp_dir/unsupported-output" -name '*.png' -print -quit | grep -q .; then
    fail "unsupported basis produced an image"
fi

JMOL_JAR="$tmp_dir/missing.jar" JMOL_TEMPLATE="$template" \
    "$renderer" "$tmp_dir/supported.molden" "$tmp_dir/missing-output" \
    > "$tmp_dir/missing.stdout" 2> "$tmp_dir/missing.stderr"
grep -q 'Jmol jar is unavailable' "$tmp_dir/missing.stderr" ||
    fail "missing Jmol warning was not emitted"

echo "PASS: render_molden_orbitals"
