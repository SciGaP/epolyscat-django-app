#!/bin/sh

set -u

executable=${EPOLYSCAT_EXECUTABLE:-/home1/06170/ampuser/apps/ePolyScat/epolyscat/bin/expanseintelmkl/ePolyScat.exe}
output_file=${EPOLYSCAT_OUTPUT_FILE:-ePolyScat.out}
launcher=${EPOLYSCAT_LAUNCHER-}
utility_bin_dir=${EPOLYSCAT_UTILITY_BIN_DIR:-/home1/06170/ampuser/apps/ePolyScat/epolyscat/bin/expanseintelmkl}
created_input=0
created_fort10=0

fail() {
    echo "Frontera ePolyScat wrapper: $1" >&2
    exit "${2:-64}"
}

cleanup() {
    if [ "$created_input" -eq 1 ]; then
        rm -f ./inp.dat
    fi
    if [ "$created_fort10" -eq 1 ]; then
        rm -f ./fort.10
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

utility_executable_name() {
    case "$1" in
    CnvMath) printf '%s\n' CnvMath.exe ;;
    CnvMatLab) printf '%s\n' CnvMatLab.exe ;;
    CnvLinFull) printf '%s\n' CnvLinFull.exe ;;
    MoldenMerge) printf '%s\n' MoldenMerge.exe ;;
    NRFPAD) printf '%s\n' NRFPAD.exe ;;
    Cube2igor) printf '%s\n' Cube2igor.exe ;;
    *) fail "unsupported utility on Frontera: ${1:-<missing>}" ;;
    esac
}

run_utility() {
    utility_name=$1
    shift
    [ "$#" -ge 1 ] || fail "$utility_name requires a utility control file and scientific data file(s)"

    control_file=$1
    shift
    minimum_data_files=1
    [ "$utility_name" != "MoldenMerge" ] || minimum_data_files=2
    [ "$#" -ge "$minimum_data_files" ] || \
        fail "$utility_name requires at least $minimum_data_files scientific data file(s)"

    [ -f "$control_file" ] && [ -r "$control_file" ] || \
        fail "utility control file is missing or unreadable: $control_file" 66
    for data_file in "$@"; do
        [ -f "$data_file" ] && [ -r "$data_file" ] || \
            fail "utility data file is missing or unreadable: $data_file" 66
    done

    utility_executable="$utility_bin_dir/$(utility_executable_name "$utility_name")"
    [ -x "$utility_executable" ] || \
        fail "utility executable is missing or not executable: $utility_executable" 69
    utility_output=${EPOLYSCAT_UTILITY_OUTPUT_FILE:-${utility_name}.out}

    trap cleanup EXIT HUP INT TERM
    if [ "$utility_name" = "Cube2igor" ]; then
        [ ! -e ./fort.10 ] || fail "refusing to overwrite an existing fort.10" 73
        ln -s "$1" ./fort.10 || fail "failed to stage Cube2igor input as fort.10" 73
        created_fort10=1
    fi

    if "$utility_executable" < "$control_file" > "$utility_output" 2>&1; then
        utility_status=0
    else
        utility_status=$?
    fi
    cat "$utility_output"
    [ "$utility_status" -eq 0 ] || \
        fail "$utility_name exited with status $utility_status; inspect $utility_output" "$utility_status"
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
UTILITY)
    utility_name=${2-}
    [ -n "$utility_name" ] || fail "UTILITY requires an application selector"
    shift 2
    run_utility "$utility_name" "$@"
    exit 0
    ;;
WORKFLOW)
    workflow_stage=${2-}
    [ "$workflow_stage" = "Analysis" ] || \
        fail "unsupported workflow stage on Frontera: ${workflow_stage:-<missing>}"
    utility_name=${3-}
    [ -n "$utility_name" ] || fail "WORKFLOW Analysis requires a utility selector"
    shift 3
    run_utility "$utility_name" "$@"
    exit 0
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
