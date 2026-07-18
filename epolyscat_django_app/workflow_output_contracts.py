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

ANALYSIS_APPLICATION_PRIORITY = (
    "CnvLinFull",
    "CnvMath",
    "CnvMatLab",
    "NRFPAD",
    "Cube2igor",
    "MoldenMerge",
)
DEFAULT_ANALYSIS_APPLICATION = ANALYSIS_APPLICATION_PRIORITY[0]


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


def _staged_input_reference(file_data):
    filename = os.path.basename(_filename(file_data))
    return f"$pt/{filename}" if filename else ""


def _serialized_file(file_data):
    if isinstance(file_data, dict):
        descriptor = file_data.get("descriptor")
        if isinstance(descriptor, dict):
            return descriptor
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
    if file_type == "plotdata":
        roles.append("plot_data_1d")
    elif file_type == "plotdata2d":
        roles.append("plot_data_2d")
    elif file_type == "mfdcs":
        roles.extend(("cross_section", "plot_data_2d"))
    elif file_type == "mfdcsfull":
        roles.append("cross_section")
    elif file_type == "mfdcsgeom":
        roles.append("cross_section_geometry")
    elif any(
        token in name
        for token in ("cross-section", "cross_section", "crosssection", "mfdcs")
    ):
        roles.append("cross_section")
    if file_type in {"matrixelements", "vibaveidy"}:
        roles.append("matrix_elements")
    if file_type in {"awavefun", "swavefun"}:
        roles.append("wavefunction")
    if file_type in {
        "dumporbbasis",
        "vieworb",
        "vieworbflux",
        "vieworbgeom",
    }:
        roles.append("orbital_grid")
    if file_type in {
        "orientasymdata",
        "orientasymeig",
        "orientasymgeom",
        "orientdata",
        "orientgeom",
    }:
        roles.append("orientation_data")
    return roles


def _files_for_role(output_files, role):
    return [
        file_data
        for file_data in output_files or []
        if role
        in (
            file_data.get("roles", [])
            if isinstance(file_data, dict) and "roles" in file_data
            else classify_output_file(file_data)
        )
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
    output_files=None,
    output_manifest=None,
    source_application,
    target_stage,
    target_application="",
    required_file_name="",
):
    """Choose the source output that satisfies one target workflow input."""

    files = list(output_manifest if output_manifest is not None else output_files or [])
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
            return {
                "status": "ready",
                "input_file_name": input_file_name,
                "expected_role": expected_role,
                "selected": _serialized_file(selected),
                "candidates": [_serialized_file(file_data) for file_data in candidates],
                "data_entry_values": {
                    "convertSource": _staged_input_reference(selected),
                    "convertFormat": (
                        "gaussian" if expected_role == "gaussian_output" else "molden"
                    ),
                },
            }
        return _missing_result(input_file_name, expected_roles[0])

    if target_stage == "Analysis" and target_application in ANALYSIS_INPUT_CONTRACTS:
        contract_input_name, expected_role = ANALYSIS_INPUT_CONTRACTS[
            target_application
        ]
        if not required_file_name or required_file_name == contract_input_name:
            input_file_name = required_file_name or contract_input_name
            candidates = _files_for_role(files, expected_role)
            if not candidates:
                return _missing_result(input_file_name, expected_role)
            return {
                "status": "ready",
                "input_file_name": input_file_name,
                "expected_role": expected_role,
                "selected": _serialized_file(candidates[0]),
                "candidates": [
                    _serialized_file(file_data) for file_data in candidates
                ],
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


def resolve_next_stage_preview(
    *, output_manifest, source_application, next_stage
):
    """Resolve the default target and inherited file for a continuation."""

    if next_stage == "ePolyScat_Run":
        binding = resolve_workflow_output_binding(
            output_manifest=output_manifest,
            source_application=source_application,
            target_stage=next_stage,
        )
        return {
            **binding,
            "next_stage": next_stage,
            "target_application": "ePolyScat",
        }

    if next_stage == "Analysis":
        fallback = None
        for application in ANALYSIS_APPLICATION_PRIORITY:
            input_file_name, _expected_role = ANALYSIS_INPUT_CONTRACTS[application]
            binding = resolve_workflow_output_binding(
                output_manifest=output_manifest,
                source_application=source_application,
                target_stage=next_stage,
                target_application=application,
                required_file_name=input_file_name,
            )
            preview = {
                **binding,
                "next_stage": next_stage,
                "target_application": application,
            }
            if fallback is None:
                fallback = preview
            if binding["status"] == "ready":
                return preview
        return fallback

    return {
        "status": "missing",
        "next_stage": next_stage,
        "target_application": "",
        "input_file_name": "",
        "expected_role": "",
        "selected": None,
        "candidates": [],
        "data_entry_values": None,
    }
