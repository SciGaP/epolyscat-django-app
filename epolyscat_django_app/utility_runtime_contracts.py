"""Shared launch contracts for the ePolyScat post-processing utilities."""


UTILITY_CONTROL_INPUT_NAMES = (
    "ePolyscat_Input_File",
    "ePolyScat_Input_File",
)

GENERIC_UTILITY_DATA_INPUT_NAMES = (
    "ePolyscat_Input_Data",
    "ePolyScat_Input_Data",
)

UTILITY_CONTRACTS = {
    "CnvMath": {
        "executable_name": "CnvMath.exe",
        "control_input_names": UTILITY_CONTROL_INPUT_NAMES,
        "data_input_names": ("BendOrient_Output",),
        "minimum_data_file_count": 1,
    },
    "CnvMatLab": {
        "executable_name": "CnvMatLab.exe",
        "control_input_names": UTILITY_CONTROL_INPUT_NAMES,
        "data_input_names": ("BendOrient_Output",),
        "minimum_data_file_count": 1,
    },
    "CnvLinFull": {
        "executable_name": "CnvLinFull.exe",
        "control_input_names": UTILITY_CONTROL_INPUT_NAMES,
        "data_input_names": ("DumpOut",),
        "minimum_data_file_count": 1,
    },
    "MoldenMerge": {
        "executable_name": "MoldenMerge.exe",
        "control_input_names": UTILITY_CONTROL_INPUT_NAMES,
        "data_input_names": ("molden.dat",),
        "minimum_data_file_count": 2,
    },
    "NRFPAD": {
        "executable_name": "NRFPAD.exe",
        "control_input_names": UTILITY_CONTROL_INPUT_NAMES,
        "data_input_names": ("Cross_Section_Input_File",),
        "minimum_data_file_count": 1,
    },
    "Cube2igor": {
        "executable_name": "Cube2igor.exe",
        "control_input_names": UTILITY_CONTROL_INPUT_NAMES,
        "data_input_names": ("Cube_Output",),
        "minimum_data_file_count": 1,
    },
}


def get_utility_contract(utility_id):
    try:
        return UTILITY_CONTRACTS[utility_id]
    except KeyError as error:
        raise ValueError(f"Unsupported ePolyScat utility: {utility_id}") from error


def _input_file_uris(input_values, input_names):
    uris = []
    for input_name in input_names:
        value = input_values.get(input_name)
        values = value if isinstance(value, (list, tuple)) else str(value or "").split(",")
        for uri in values:
            uri = str(uri).strip()
            if uri and uri not in uris:
                uris.append(uri)
    return uris


def _first_populated_input_name(input_values, input_names):
    return next((name for name in input_names if input_values.get(name)), None)


def resolve_utility_argument_input_names(utility_id, input_values):
    contract = get_utility_contract(utility_id)
    input_names = set()

    control_input_name = _first_populated_input_name(
        input_values,
        contract["control_input_names"],
    )
    if control_input_name:
        input_names.add(control_input_name)

    specific_names = tuple(
        name for name in contract["data_input_names"] if input_values.get(name)
    )
    specific_uris = _input_file_uris(input_values, specific_names)
    input_names.update(specific_names)

    for generic_name in GENERIC_UTILITY_DATA_INPUT_NAMES:
        if not input_values.get(generic_name):
            continue
        generic_uris = _input_file_uris(input_values, (generic_name,))
        if not specific_names or any(uri not in specific_uris for uri in generic_uris):
            input_names.add(generic_name)
            break

    return input_names


def resolve_utility_argument_input_order(utility_id, input_values):
    contract = get_utility_contract(utility_id)
    selected_names = resolve_utility_argument_input_names(
        utility_id,
        input_values,
    )
    ordered_names = []

    for input_name in contract["control_input_names"]:
        if input_name in selected_names:
            ordered_names.append(input_name)
            break

    for input_name in contract["data_input_names"]:
        if input_name in selected_names:
            ordered_names.append(input_name)

    for input_name in GENERIC_UTILITY_DATA_INPUT_NAMES:
        if input_name in selected_names:
            ordered_names.append(input_name)

    return tuple(ordered_names)


def validate_utility_input_values(utility_id, input_values):
    contract = get_utility_contract(utility_id)
    argument_input_names = resolve_utility_argument_input_names(
        utility_id,
        input_values,
    )
    control_input_names = set(contract["control_input_names"])
    if not argument_input_names.intersection(control_input_names):
        raise ValueError(f"{utility_id} requires a utility control file")

    data_input_names = argument_input_names.difference(control_input_names)
    data_uris = _input_file_uris(input_values, data_input_names)
    minimum_file_count = contract["minimum_data_file_count"]
    if len(data_uris) < minimum_file_count:
        preferred_name = contract["data_input_names"][0]
        raise ValueError(
            f"{utility_id} requires {preferred_name} with at least "
            f"{minimum_file_count} scientific data file(s)"
        )
