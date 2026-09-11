#!/bin/bash

set -uo pipefail

warn() {
    echo "Molden orbital renderer: $1" >&2
}

fail() {
    warn "$1"
    exit "${2:-1}"
}

molden_file=${1-}
output_dir=${2:-$PWD}

[[ -n "$molden_file" ]] || fail "usage: $0 molden_file [output_dir]" 64
[[ -r "$molden_file" ]] || fail "input file is missing or unreadable: $molden_file" 66
mkdir -p "$output_dir" || fail "cannot create output directory: $output_dir" 73

molden_dir=$(cd -- "$(dirname -- "$molden_file")" && pwd) ||
    fail "cannot resolve input directory: $molden_file" 66
molden_file="$molden_dir/$(basename -- "$molden_file")"
output_dir=$(cd -- "$output_dir" && pwd) ||
    fail "cannot resolve output directory: $output_dir" 73

unsupported_marker=$(awk '
    {
        line = $0
        sub(/\r$/, "", line)
        if (line ~ /^\[GTO\]/ || line ~ /^\[GTO\][[:space:]]/) {
            in_gto = 1
        } else if (line ~ /^\[/) {
            in_gto = 0
        }
        if (line ~ /^\[([0-9]+[A-Za-z])+\]$/ && line ~ /[GHIghi]/) {
            print line
            exit
        }
        if (in_gto && line ~ /^[[:space:]]*[GHIghi][[:space:]]+[0-9]+([[:space:]]|$)/) {
            sub(/^[[:space:]]*/, "", line)
            print line
            exit
        }
    }
' "$molden_file")

if [[ -n "$unsupported_marker" ]]; then
    warn "orbital previews skipped: unsupported angular momentum marker $unsupported_marker"
    exit 0
fi

script_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
app_dir=$(cd -- "$script_dir/.." && pwd)
jmol_jar=${JMOL_JAR:-/home1/06170/ampuser/AMP_apps/Jmol-SwingJS/dist/JmolDataD.jar}
jmol_template=${JMOL_TEMPLATE:-$app_dir/templates/print_frontier_orbitals.template}
java_cmd=${JMOL_JAVA:-java}

if [[ ! -r "$jmol_jar" ]]; then
    warn "orbital previews skipped: Jmol jar is unavailable: $jmol_jar"
    exit 0
fi
if [[ ! -r "$jmol_template" ]]; then
    warn "orbital previews skipped: Jmol template is unavailable: $jmol_template"
    exit 0
fi
if ! command -v "$java_cmd" >/dev/null 2>&1; then
    warn "orbital previews skipped: Java is unavailable: $java_cmd"
    exit 0
fi

escape_sed_replacement() {
    printf '%s' "$1" | sed -e 's/[\\&|]/\\&/g'
}

escaped_molden=$(escape_sed_replacement "$molden_file")
escaped_output=$(escape_sed_replacement "$output_dir")
work_dir=$(mktemp -d "$output_dir/.jmol-orbitals.XXXXXX") ||
    fail "cannot create temporary Jmol directory" 73
trap 'rm -rf "$work_dir"' EXIT

export DISPLAY=""
views=(default x y z)
for view in "${views[@]}"; do
    if [[ "$view" == "default" ]]; then
        rotation=""
        suffix=""
    else
        rotation="rotate $view 90;"
        suffix="_$view"
    fi

    escaped_rotation=$(escape_sed_replacement "$rotation")
    escaped_suffix=$(escape_sed_replacement "$suffix")
    jmol_script="$work_dir/frontier_orbitals_$view.jmol"
    jmol_log="$work_dir/frontier_orbitals_$view.log"

    sed \
        -e "s|@molden_file|$escaped_molden|g" \
        -e "s|@output_dir|$escaped_output|g" \
        -e "s|@rotate|$escaped_rotation|g" \
        -e "s|@dir|$escaped_suffix|g" \
        "$jmol_template" > "$jmol_script" ||
        fail "could not create Jmol script for view $view" 74

    if ! "$java_cmd" -jar "$jmol_jar" -g800x600 -ions "$jmol_script" > "$jmol_log" 2>&1; then
        cp "$jmol_script" "$output_dir/jmol_orbitals_failed_$view.jmol"
        cp "$jmol_log" "$output_dir/jmol_orbitals_failed_$view.log"
        fail "Jmol failed for view $view; diagnostic files were preserved" 70
    fi

    for orbital in homo lumo; do
        image="$output_dir/orbital_${orbital}${suffix}.png"
        if [[ ! -s "$image" ]]; then
            cp "$jmol_script" "$output_dir/jmol_orbitals_failed_$view.jmol"
            cp "$jmol_log" "$output_dir/jmol_orbitals_failed_$view.log"
            fail "Jmol did not create $image; diagnostic files were preserved" 70
        fi
    done
done

echo "Molden orbital renderer: generated HOMO and LUMO previews in four orientations"
