#!/bin/sh

set -u

bin_dir=${EPOLYSCAT_BIN_DIR:-/home/ampuser/apps/ePolyScat/epolyscat/bin}
runtime_dir=${EPOLYSCAT_RUNTIME_DIR:-$bin_dir/expanseintelmkl}
executable=${EPOLYSCAT_EXECUTABLE:-$runtime_dir/ePolyScat.exe}
utility_bin_dir=${EPOLYSCAT_UTILITY_BIN_DIR:-$runtime_dir}
output_file=${EPOLYSCAT_OUTPUT_FILE:-ePolyScat.out}
launcher=${EPOLYSCAT_LAUNCHER-}
created_input=0
created_fort10=0
created_fort21=0
created_fort22=0

fail() {
    echo "Expanse ePolyScat controller: $1" >&2
    exit "${2:-64}"
}

cleanup() {
    if [ "$created_input" -eq 1 ]; then
        rm -f ./inp.dat
    fi
    if [ "$created_fort10" -eq 1 ]; then
        rm -f ./fort.10
    fi
    if [ "$created_fort21" -eq 1 ]; then
        rm -f ./fort.21
    fi
    if [ "$created_fort22" -eq 1 ]; then
        rm -f ./fort.22
    fi
}

command_exists() {
    command -v "$1" >/dev/null 2>&1
}

setup_expanse_runtime() {
    if [ "${EPOLYSCAT_SKIP_MODULE_SETUP:-0}" = "1" ]; then
        return
    fi

    if ! command_exists module; then
        if [ -r /etc/profile.d/modules.sh ]; then
            . /etc/profile.d/modules.sh
        elif [ -r /etc/profile.d/lmod.sh ]; then
            . /etc/profile.d/lmod.sh
        fi
    fi
    command_exists module || fail "the Expanse module command is unavailable" 69

    module --force purge >/dev/null 2>&1 || fail "failed to purge the module environment" 69
    module load \
        cpu/0.15.4 \
        intel/19.1.1.217 \
        intel-mpi/2019.8.254 \
        intel-mkl/2019.1.144 \
        slurm >/dev/null 2>&1 || \
        fail "failed to load the validated Expanse Intel runtime" 69
}

escape_sed_replacement() {
    printf '%s' "$1" | sed 's/[\\&|]/\\&/g'
}

stage_input() {
    source_file=$1

    if ! LC_ALL=C grep -Eq '\$(pt|po|pd)/' "$source_file"; then
        cp "$source_file" ./inp.dat
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

slurm_task_count() {
    task_count=${SLURM_NTASKS:-${SLURM_NPROCS:-1}}
    case "$task_count" in
    ""|0|*[!0-9]*)
        fail "invalid SLURM task count: ${task_count:-<missing>}"
        ;;
    esac
    printf '%s\n' "$task_count"
}

run_epolyscat_executable() {
    task_count=$(slurm_task_count)

    if [ -n "$launcher" ]; then
        command_exists "$launcher" || fail "launcher is unavailable: $launcher" 69
        "$launcher" -bootstrap slurm -np "$task_count" "$executable"
    elif [ -n "${SLURM_JOB_ID-}" ]; then
        command_exists mpirun || fail "Expanse launcher is unavailable: mpirun" 69
        mpirun -bootstrap slurm -np "$task_count" "$executable"
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
    *) fail "unsupported utility on Expanse: ${1:-<missing>}" ;;
    esac
}

run_utility() {
    utility_name=$1
    shift
    [ "$#" -ge 1 ] || \
        fail "$utility_name requires a utility control file and scientific data file(s)"

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
    if [ "$utility_name" = "NRFPAD" ] && [ "$#" -gt 2 ]; then
        fail "NRFPAD accepts one FLNNuLP file and one optional HNu file"
    fi

    setup_expanse_runtime
    utility_executable="$utility_bin_dir/$(utility_executable_name "$utility_name")"
    [ -x "$utility_executable" ] || \
        fail "utility executable is missing or not executable: $utility_executable" 69
    utility_output=${EPOLYSCAT_UTILITY_OUTPUT_FILE:-${utility_name}.out}

    trap cleanup EXIT HUP INT TERM
    if [ "$utility_name" = "Cube2igor" ]; then
        case "$1" in
        fort.10|./fort.10)
            ;;
        *)
            [ ! -e ./fort.10 ] || fail "refusing to overwrite an existing fort.10" 73
            ln -s "$1" ./fort.10 || fail "failed to stage Cube2igor input as fort.10" 73
            created_fort10=1
            ;;
        esac
    fi
    if [ "$utility_name" = "NRFPAD" ]; then
        nrfpad_fln_file=$1
        case "$nrfpad_fln_file" in
        fort.21|./fort.21)
            ;;
        *)
            [ ! -e ./fort.21 ] || fail "refusing to overwrite an existing fort.21" 73
            ln -s "$nrfpad_fln_file" ./fort.21 || \
                fail "failed to stage NRFPAD FLNNuLP input as fort.21" 73
            created_fort21=1
            ;;
        esac
        if [ "$#" -eq 2 ]; then
            nrfpad_hnu_file=$2
            case "$nrfpad_hnu_file" in
            fort.22|./fort.22)
                ;;
            *)
                [ ! -e ./fort.22 ] || fail "refusing to overwrite an existing fort.22" 73
                ln -s "$nrfpad_hnu_file" ./fort.22 || \
                    fail "failed to stage NRFPAD HNu input as fort.22" 73
                created_fort22=1
                ;;
            esac
        fi
    fi

    if "$utility_executable" < "$control_file" > "$utility_output" 2>&1; then
        utility_status=0
    else
        utility_status=$?
    fi
    cat "$utility_output"
    [ "$utility_status" -eq 0 ] || \
        fail "$utility_name exited with status $utility_status; inspect $utility_output" "$utility_status"
    if LC_ALL=C grep -Eiq \
        'ABSTOP|Abnormal Ending|forrtl: severe|Segmentation fault|Floating point exception' \
        "$utility_output"; then
        fail "$utility_name reported a scientific runtime failure; inspect $utility_output" 70
    fi
}

run_delegate() {
    delegate_name=$1
    shift
    if "$@"; then
        return
    else
        delegate_status=$?
    fi
    fail "$delegate_name exited with status $delegate_status" "$delegate_status"
}

run_gaussian() {
    [ "$#" -ge 1 ] || fail "Gaussian16 requires a staged input filename"
    gaussian_input=$1
    gaussian_runner="$bin_dir/eps_g16.sh"
    [ -x "$gaussian_runner" ] || \
        fail "Gaussian16 runner is missing or not executable: $gaussian_runner" 69
    run_delegate Gaussian16 "$gaussian_runner" "$@"

    gaussian_input_name=${gaussian_input##*/}
    gaussian_output_name=${gaussian_input_name%.*}.log
    [ -s "$gaussian_output_name" ] || \
        fail "Gaussian16 did not produce a non-empty output log: $gaussian_output_name" 70
    LC_ALL=C grep -q "Normal termination of Gaussian" "$gaussian_output_name" || \
        fail "Gaussian16 output does not contain a normal termination marker: $gaussian_output_name" 70
}

run_openmolcas() {
    [ "$#" -eq 1 ] || fail "OpenMolcas requires exactly one staged input filename"
    openmolcas_runner="$bin_dir/eps_openmolcas.sh"
    [ -x "$openmolcas_runner" ] || \
        fail "OpenMolcas runner is missing or not executable: $openmolcas_runner" 69
    run_delegate OpenMolcas "$openmolcas_runner" "$1" -i=yes
}

run_epolyscat() {
    input_file=$1

    [ -f "$input_file" ] && [ -r "$input_file" ] || \
        fail "input file is missing or unreadable: $input_file" 66
    [ "$output_file" != "inp.dat" ] && [ "$output_file" != "./inp.dat" ] || \
        fail "output file must not be inp.dat" 73

    setup_expanse_runtime
    [ -x "$executable" ] || \
        fail "scientific executable is missing or not executable: $executable" 69

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

    if run_epolyscat_executable > "$output_file" 2>&1; then
        run_status=0
    else
        run_status=$?
    fi
    cat "$output_file"

    [ "$run_status" -eq 0 ] || \
        fail "launcher or scientific executable exited with status $run_status; inspect $output_file" "$run_status"
    if LC_ALL=C grep -q "Abnormal Ending" "$output_file"; then
        fail "calculation reported Abnormal Ending; inspect $output_file" 70
    fi
    LC_ALL=C grep -q "Finalize" "$output_file" || \
        fail "calculation did not reach Finalize; inspect $output_file" 70
}

run_type=${1-}

case "$run_type" in
MODULE)
    module_name=${2-}
    case "$module_name" in
    ePolyScat)
        [ "$#" -le 3 ] || fail "unexpected arguments after MODULE ePolyScat"
        input_file=${3-${EPOLYSCAT_INPUT_FILE-}}
        [ -n "$input_file" ] || input_file=ep_input
        run_epolyscat "$input_file"
        ;;
    Gaussian16)
        shift 2
        run_gaussian "$@"
        ;;
    OpenMolcas)
        shift 2
        run_openmolcas "$@"
        ;;
    *)
        fail "unsupported module on Expanse: ${module_name:-<missing>}"
        ;;
    esac
    ;;
UTILITY)
    utility_name=${2-}
    [ -n "$utility_name" ] || fail "UTILITY requires an application selector"
    shift 2
    run_utility "$utility_name" "$@"
    ;;
WORKFLOW)
    workflow_stage=${2-}
    case "$workflow_stage" in
    Data_Gen)
        workflow_application=${3-}
        [ -n "$workflow_application" ] || \
            fail "WORKFLOW Data_Gen requires an application selector"
        shift 3
        case "$workflow_application" in
        Gaussian16) run_gaussian "$@" ;;
        OpenMolcas) run_openmolcas "$@" ;;
        *) fail "unsupported Data_Gen application on Expanse: $workflow_application" ;;
        esac
        ;;
    ePolyScat_Run)
        [ "$#" -le 3 ] || fail "unexpected arguments after WORKFLOW ePolyScat_Run"
        input_file=${3-${EPOLYSCAT_INPUT_FILE-}}
        [ -n "$input_file" ] || input_file=ep_input
        run_epolyscat "$input_file"
        ;;
    Analysis)
        utility_name=${3-}
        [ -n "$utility_name" ] || \
            fail "WORKFLOW Analysis requires a utility selector"
        shift 3
        run_utility "$utility_name" "$@"
        ;;
    *)
        fail "unsupported workflow stage on Expanse: ${workflow_stage:-<missing>}"
        ;;
    esac
    ;;
"")
    fail "run type is required"
    ;;
*)
    fail "unsupported run type: $run_type"
    ;;
esac

exit 0
