"""Scientific output contracts used to connect ePolyScat workflow stages."""

import os
import re


ANALYSIS_INPUT_CONTRACTS = {
    "CnvMath": ("BendOrient_Output", "bend_orient"),
    "CnvMatLab": ("BendOrient_Output", "bend_orient"),
    "CnvLinFull": ("DumpOut", "dump_out"),
    "MoldenMerge": ("molden.dat", "molden"),
    "NRFPAD": ("Cross_Section_Input_File", "orient_ncro"),
    "Cube2igor": ("Cube_Output", "cube"),
}


def _value(file_data, *names):
    for name in names:
        if isinstance(file_data, dict):
            value = file_data.get(name)
        else:
            value = getattr(file_data, name, None)
        if value not in (None, ""):
            return value
    return ""


def _normalized_token(value):
    return re.sub(r"[^a-z0-9]", "", str(value or "").lower())


def _filename(file_data):
    return str(_value(file_data, "name", "filename", "fileName") or "")


def _serialized_file(file_data):
    if isinstance(file_data, dict):
        return file_data
    try:
        return {
            key: value
            for key, value in vars(file_data).items()
            if not key.startswith("_")
        }
    except TypeError:
        return {"name": _filename(file_data)}


def classify_output_file(file_data):
    """Return semantic roles for an Airavata output descriptor."""

    name = os.path.basename(_filename(file_data)).lower()
    file_type = _normalized_token(
        _value(file_data, "fileType", "file_type", "type", "outputType")
    )
    if file_type in {"stdout", "stderr"} or name.endswith(
        (".stdout", ".stderr", ".slurm")
    ):
        return []

    roles = []
    if file_type in {"gaussianlog", "gaussianoutput"} or (
        name == "gaussian.log"
        or name == "fort.7"
        or re.search(r"\.(g03|g09|g16)$", name)
        or ("gaussian" in name and name.endswith((".log", ".out")))
    ):
        roles.append("gaussian_output")
    if file_type in {"molden", "moldenoutput", "moldenfile"} or (
        "molden" in name or name.endswith(".molden")
    ):
        roles.append("molden")
    if file_type in {"dumpout", "dumpidy", "dumpidyall"} or (
        "dumpout" in name or "dumpidy" in name
    ):
        roles.append("dump_out")
    if file_type in {"bendorient", "bendorientoutput"} or "bendorient" in name:
        roles.append("bend_orient")
    if file_type in {"orientncro", "orientncrooutput", "orinetncro"} or (
        "orientncro" in name or "orinetncro" in name
    ):
        roles.append("orient_ncro")
    if file_type in {"cube", "cubeoutput", "gaussiancube"} or name.endswith(
        (".cube", ".cub")
    ):
        roles.append("cube")
    return roles


def _files_for_role(output_files, role):
    return [
        file_data
        for file_data in output_files or []
        if role in classify_output_file(file_data)
    ]


def _missing_result(input_file_name, expected_role):
    return {
        "status": "missing",
        "input_file_name": input_file_name,
        "expected_role": expected_role,
        "selected": None,
        "candidates": [],
        "data_entry_values": None,
    }

def resolve_workflow_output_binding(
    *,
    output_files,
    source_application,
    target_stage,
    target_application="",
    required_file_name="",
):
    """Choose the source output that satisfies one target workflow input."""

    files = list(output_files or [])
    source_application = str(source_application or "")
    target_stage = str(target_stage or "")

    if target_stage == "ePolyScat_Run":
        expected_roles = (
            ("gaussian_output",)
            if source_application == "Gaussian16"
            else ("molden",)
            if source_application == "OpenMolcas"
            else ("molden", "gaussian_output")
        )
        input_file_name = "ePolyScat_Input_Data"
        for expected_role in expected_roles:
            candidates = _files_for_role(files, expected_role)
            if not candidates:
                continue
            selected = candidates[0]
            filename = _filename(selected)
            return {
                "status": "ready",
                "input_file_name": input_file_name,
                "expected_role": expected_role,
                "selected": _serialized_file(selected),
                "candidates": [_serialized_file(file_data) for file_data in candidates],
                "data_entry_values": {
                    "convertSource": filename,
                    "convertFormat": (
                        "gaussian" if expected_role == "gaussian_output" else "molden"
                    ),
                },
            }
        return _missing_result(input_file_name, expected_roles[0])

    if target_stage == "Analysis" and target_application in ANALYSIS_INPUT_CONTRACTS:
        input_file_name, expected_role = ANALYSIS_INPUT_CONTRACTS[target_application]
        input_file_name = required_file_name or input_file_name
        candidates = _files_for_role(files, expected_role)
        if not candidates:
            return _missing_result(input_file_name, expected_role)
        return {
            "status": "ready",
            "input_file_name": input_file_name,
            "expected_role": expected_role,
            "selected": _serialized_file(candidates[0]),
            "candidates": [_serialized_file(file_data) for file_data in candidates],
            "data_entry_values": None,
        }

    required_name = str(required_file_name or "").lower()
    candidates = [
        file_data for file_data in files if _filename(file_data).lower() == required_name
    ]
    if not candidates:
        return _missing_result(required_file_name, "exact_filename")
    return {
        "status": "ready",
        "input_file_name": required_file_name,
        "expected_role": "exact_filename",
        "selected": _serialized_file(candidates[0]),
        "candidates": [_serialized_file(file_data) for file_data in candidates],
        "data_entry_values": None,
    }
