#!/bin/sh

set -u

executable=${EPOLYSCAT_EXECUTABLE:-/home1/06170/ampuser/apps/ePolyScat/epolyscat/bin/expanseintelmkl/ePolyScat.exe}
output_file=${EPOLYSCAT_OUTPUT_FILE:-ePolyScat.out}
launcher=${EPOLYSCAT_LAUNCHER-}
created_input=0

fail() {
    echo "Frontera ePolyScat wrapper: $1" >&2
    exit "${2:-64}"
}

cleanup() {
    if [ "$created_input" -eq 1 ]; then
        rm -f ./inp.dat
    fi
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

escape_sed_replacement() {
    printf '%s' "$1" | sed 's/[\\&|]/\\&/g'
}

stage_input() {
    source_file=$1

    if ! LC_ALL=C grep -Eq '\$(pt|po|pd)/' "$source_file"; then
        cp -- "$source_file" ./inp.dat
        return
    fi

    pt_path=${pt:-.}
    po_path=${po:-.}
    pd_path=${pd:-.}
    pt_path=${pt_path%/}
    po_path=${po_path%/}
    pd_path=${pd_path%/}

    pt_replacement=$(escape_sed_replacement "$pt_path")
    po_replacement=$(escape_sed_replacement "$po_path")
    pd_replacement=$(escape_sed_replacement "$pd_path")

    LC_ALL=C sed \
        -e "s|\\\$pt/|$pt_replacement/|g" \
        -e "s|\\\$po/|$po_replacement/|g" \
        -e "s|\\\$pd/|$pd_replacement/|g" \
        "$source_file" > ./inp.dat
}

run_epolyscat() {
    if [ -n "$launcher" ]; then
        command_exists "$launcher" || fail "launcher is unavailable: $launcher" 69
        "$launcher" "$executable"
    elif [ -n "${SLURM_JOB_ID-}" ]; then
        command_exists ibrun || fail "Frontera launcher is unavailable: ibrun" 69
        ibrun "$executable"
    else
        "$executable"
    fi
}

run_type=${1-}

case "$run_type" in
MODULE)
    module_name=${2-}
    [ "$module_name" = "ePolyScat" ] || fail "unsupported module on Frontera: ${module_name:-<missing>}"
    [ "$#" -le 3 ] || fail "unexpected arguments after MODULE ePolyScat"
    input_file=${3-${EPOLYSCAT_INPUT_FILE-}}
    [ -n "$input_file" ] || input_file=ep_input
    ;;
UTILITY|WORKFLOW)
    fail "unsupported run type on the Frontera ePolyScat deployment: $run_type"
    ;;
"")
    fail "an input file or MODULE ePolyScat selector is required"
    ;;
*)
    [ "$#" -eq 1 ] || fail "direct invocation accepts exactly one input file"
    input_file=$1
    ;;
esac

[ -f "$input_file" ] && [ -r "$input_file" ] || \
    fail "input file is missing or unreadable: $input_file" 66
[ -x "$executable" ] || fail "scientific executable is missing or not executable: $executable" 69
[ "$output_file" != "inp.dat" ] && [ "$output_file" != "./inp.dat" ] || \
    fail "output file must not be inp.dat" 73

case "$input_file" in
inp.dat|./inp.dat)
    ;;
*)
    [ ! -e ./inp.dat ] || fail "refusing to overwrite an existing inp.dat" 73
    stage_input "$input_file" || fail "failed to stage input as inp.dat" 73
    created_input=1
    ;;
esac

trap cleanup EXIT HUP INT TERM

export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export MPI_MEMMAP_OFF=${MPI_MEMMAP_OFF:-1}
export MPI_DSM_DISTRIBUTE=${MPI_DSM_DISTRIBUTE:-1}

if run_epolyscat > "$output_file" 2>&1; then
    run_status=0
else
    run_status=$?
fi

cat "$output_file"

[ "$run_status" -eq 0 ] || \
    fail "launcher or scientific executable exited with status $run_status; inspect $output_file" "$run_status"

if grep -q "Abnormal Ending" "$output_file"; then
    fail "calculation reported Abnormal Ending; inspect $output_file" 70
fi

grep -q "Finalize" "$output_file" || \
    fail "calculation did not reach Finalize; inspect $output_file" 70

exit 0
