"""Presentation metadata for heterogeneous scientific output files."""

import os
import re
import shlex

from .workflow_output_contracts import classify_output_file


PLOT_CONTRACTS = {
    "plotdata": {
        "dimension": 1,
        "x_axis": "0",
        "y_axes": "1",
        "flags": "-linY",
    },
    "plotdata2d": {
        "dimension": 2,
        "x_axis": "0,1",
        "y_axes": "2",
        "flags": "",
    },
    "mfdcs": {
        "dimension": 2,
        "x_axis": "0,1",
        "y_axes": "2",
        "flags": "",
    },
}

BINARY_FILE_EXTENSIONS = {
    ".7z",
    ".bin",
    ".chk",
    ".dcd",
    ".gbw",
    ".gz",
    ".h5",
    ".hdf5",
    ".jpeg",
    ".jpg",
    ".pdf",
    ".png",
    ".rwf",
    ".tar",
    ".tif",
    ".tiff",
    ".zip",
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


def _filename(file_data):
    return os.path.basename(
        str(_value(file_data, "name", "filename", "fileName") or "")
    )


def parse_file_name_declarations(contents):
    """Return output basename-to-type bindings declared by ePolyScat FileName."""

    if isinstance(contents, bytes):
        contents = contents.decode("utf-8", errors="replace")

    declarations = {}
    for line in str(contents or "").splitlines():
        try:
            tokens = shlex.split(line, comments=True, posix=True)
        except ValueError:
            continue
        if len(tokens) < 3 or tokens[0].lower() != "filename":
            continue

        filename = os.path.basename(tokens[2].replace("\\", "/"))
        if filename:
            declarations[filename] = tokens[1]
    return declarations


def _declared_file_type(filename, declared_file_types):
    if filename in declared_file_types:
        return declared_file_types[filename]
    lower_filename = filename.lower()
    for declared_name, file_type in declared_file_types.items():
        if str(declared_name).lower() == lower_filename:
            return file_type
    return ""


def is_viewable_file(filename):
    extension = os.path.splitext(filename.lower())[1]
    return extension not in BINARY_FILE_EXTENSIONS


def annotate_output_files(output_files, declared_file_types=None):
    """Preserve output descriptors while adding display and plotting contracts."""

    declared_file_types = dict(declared_file_types or {})
    annotated = []
    for file_data in output_files or []:
        descriptor = _serialize_descriptor(file_data)
        filename = _filename(file_data)
        declared_type = _declared_file_type(filename, declared_file_types)
        file_type = declared_type or str(
            _value(file_data, "fileType", "file_type", "type", "outputType")
            or ""
        )

        classified_descriptor = dict(descriptor)
        classified_descriptor["name"] = filename
        classified_descriptor["fileType"] = file_type
        semantic_roles = list(classify_output_file(classified_descriptor))
        plot_contract = PLOT_CONTRACTS.get(_normalized_token(file_type))

        descriptor.setdefault("name", filename)
        descriptor["file_type"] = file_type
        descriptor["semantic_roles"] = semantic_roles
        descriptor["viewable"] = is_viewable_file(filename)
        descriptor["plottable"] = plot_contract is not None
        descriptor["plot_contract"] = (
            dict(plot_contract) if plot_contract is not None else None
        )
        annotated.append(descriptor)
    return annotated
