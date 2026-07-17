"""Scientific output manifests and success verification rules."""

import os
import re

from .workflow_output_contracts import (
    ANALYSIS_INPUT_CONTRACTS,
    classify_output_file,
)


EPOLYSCAT_RESULT_ROLES = {
    "bend_orient",
    "cross_section",
    "cube",
    "dump_out",
    "matrix_elements",
    "orbital_grid",
    "orientation_data",
    "orient_ncro",
    "wavefunction",
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


def _serialize_descriptor(file_data):
    if isinstance(file_data, dict):
        return dict(file_data)
    try:
        return {
            key: value
            for key, value in vars(file_data).items()
            if not key.startswith("_")
        }
    except TypeError:
        return {"name": str(file_data)}


def _normalized_token(value):
    return re.sub(r"[^a-z0-9]", "", str(value or "").lower())


def _is_scheduler_artifact(name, file_type):
    normalized_type = _normalized_token(file_type)
    lower_name = str(name or "").lower()
    return normalized_type in {"stdout", "stderr"} or lower_name.endswith(
        (".stdout", ".stderr", ".slurm")
    )


def build_output_manifest(output_files):
    """Normalize heterogeneous Airavata output descriptors without losing them."""

    manifest = []
    for file_data in output_files or []:
        name = os.path.basename(
            str(_value(file_data, "name", "filename", "fileName") or "")
        )
        file_type = str(
            _value(file_data, "fileType", "file_type", "type", "outputType")
            or ""
        )
        roles = list(classify_output_file(file_data))
        scheduler_artifact = _is_scheduler_artifact(name, file_type)
        manifest.append(
            {
                "name": name,
                "data_product_uri": str(
                    _value(
                        file_data,
                        "data_product_uri",
                        "dataProductURI",
                        "data-product-uri",
                        "productUri",
                        "product_uri",
                        "value",
                    )
                    or ""
                ),
                "file_type": file_type,
                "roles": roles,
                "scheduler_artifact": scheduler_artifact,
                "scientific": bool(roles) and not scheduler_artifact,
                "descriptor": _serialize_descriptor(file_data),
            }
        )
    return manifest


def _result(application, status, reason, message, evidence_files=None):
    return {
        "application": application,
        "status": status,
        "reason": reason,
        "message": message,
        "evidence_files": list(evidence_files or []),
    }


def _read_text(entry, read_text):
    if read_text is None:
        return None
    try:
        contents = read_text(entry)
    except Exception:
        return None
    if contents is None:
        return None
    if isinstance(contents, bytes):
        return contents.decode("utf-8", errors="replace")
    return str(contents)


def _texts(entries, read_text):
    readable = []
    for entry in entries:
        contents = _read_text(entry, read_text)
        if contents is not None:
            readable.append((entry, contents))
    return readable


def _entries_with_role(manifest, role):
    return [entry for entry in manifest if role in entry.get("roles", [])]


def _log_entries(manifest, application_token):
    entries = []
    for entry in manifest:
        name = entry.get("name", "").lower()
        if entry.get("roles") == ["molden"]:
            continue
        if (
            application_token in _normalized_token(name)
            or name.endswith((".log", ".out", ".stdout"))
            or _normalized_token(entry.get("file_type")) == "stdout"
        ):
            entries.append(entry)
    return entries


def verify_scientific_outputs(application, manifest, read_text=None):
    """Verify scientific completion independently of scheduler completion."""

    application = str(application or "")
    token = _normalized_token(application)
    files = list(manifest or [])

    if token in {"gaussian", "gaussian16"}:
        outputs = _entries_with_role(files, "gaussian_output")
        if not outputs:
            return _result(
                application,
                "incomplete",
                "missing_scientific_output",
                "Gaussian did not produce a recognized output log.",
            )
        readable = _texts(outputs, read_text)
        for entry, contents in readable:
            if "normal termination of gaussian" in contents.lower():
                return _result(
                    application,
                    "verified",
                    "",
                    "Gaussian reported normal termination.",
                    [entry["name"]],
                )
        if readable:
            return _result(
                application,
                "failed",
                "success_marker_missing",
                "Gaussian output does not contain a normal termination marker.",
                [entry["name"] for entry, _contents in readable],
            )
        return _result(
            application,
            "unverified",
            "evidence_unavailable",
            "Gaussian output exists, but its contents could not be read.",
            [entry["name"] for entry in outputs],
        )

    if token == "openmolcas":
        molden_outputs = _entries_with_role(files, "molden")
        if not molden_outputs:
            return _result(
                application,
                "incomplete",
                "missing_scientific_output",
                "OpenMolcas did not produce a Molden orbital file.",
            )
        logs = _log_entries(files, "molcas")
        readable = _texts(logs, read_text)
        for entry, contents in readable:
            normalized_contents = contents.lower()
            if (
                "_rc_all_is_well_" in normalized_contents
                or "happy landing" in normalized_contents
            ):
                return _result(
                    application,
                    "verified",
                    "",
                    "OpenMolcas completed and produced a Molden orbital file.",
                    [entry["name"], molden_outputs[0]["name"]],
                )
        if readable:
            return _result(
                application,
                "failed",
                "success_marker_missing",
                "OpenMolcas output does not contain a successful landing marker.",
                [entry["name"] for entry, _contents in readable],
            )
        return _result(
            application,
            "unverified",
            "evidence_unavailable",
            "OpenMolcas produced a Molden file, but its log could not be verified.",
            [entry["name"] for entry in molden_outputs],
        )

    if token == "epolyscat":
        scientific_results = [
            entry
            for entry in files
            if EPOLYSCAT_RESULT_ROLES.intersection(entry.get("roles", []))
        ]
        if not scientific_results:
            return _result(
                application,
                "incomplete",
                "missing_scientific_output",
                "ePolyScat did not produce a recognized scientific result file.",
            )
        logs = _log_entries(files, "epolyscat")
        readable = _texts(logs, read_text)
        for entry, contents in readable:
            if "abnormal ending" in contents.lower():
                return _result(
                    application,
                    "failed",
                    "abnormal_ending",
                    "ePolyScat reported an abnormal ending.",
                    [entry["name"]],
                )
            if "End EDCS" in contents and "Finalize" in contents:
                return _result(
                    application,
                    "verified",
                    "",
                    "ePolyScat reached EDCS completion and finalization.",
                    [entry["name"], scientific_results[0]["name"]],
                )
        if readable:
            return _result(
                application,
                "failed",
                "success_marker_missing",
                "ePolyScat output is missing EDCS completion or finalization markers.",
                [entry["name"] for entry, _contents in readable],
            )
        return _result(
            application,
            "unverified",
            "evidence_unavailable",
            "ePolyScat result files exist, but its execution log could not be read.",
            [entry["name"] for entry in scientific_results],
        )

    utility_contract = ANALYSIS_INPUT_CONTRACTS.get(application)
    if utility_contract:
        _input_name, expected_role = utility_contract
        outputs = _entries_with_role(files, expected_role)
        if outputs:
            return _result(
                application,
                "verified",
                "",
                f"{application} produced its expected scientific output.",
                [outputs[0]["name"]],
            )
        return _result(
            application,
            "incomplete",
            "missing_scientific_output",
            f"{application} did not produce its expected scientific output.",
        )

    return _result(
        application,
        "unverified",
        "unsupported_application",
        "No scientific verification rule is registered for this application.",
    )
