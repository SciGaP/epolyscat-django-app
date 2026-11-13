"""Audit the shared Airavata wrapper used by ePolyScat modules and utilities."""

import json

from epolyscat_django_app import utility_runtime_contracts

EXPECTED_UTILITIES = tuple(utility_runtime_contracts.UTILITY_CONTRACTS)

UTILITY_REQUIRED_INPUTS = {
    utility_id: contract["data_input_names"]
    for utility_id, contract in utility_runtime_contracts.UTILITY_CONTRACTS.items()
}

GENERIC_POSTPROCESSING_INPUTS = utility_runtime_contracts.GENERIC_UTILITY_DATA_INPUT_NAMES


def _value(instance, *names, default=None):
    for name in names:
        if isinstance(instance, dict):
            value = instance.get(name)
        else:
            value = getattr(instance, name, None)
        if value not in (None, ""):
            return value
    return default


def _metadata_values(metadata):
    if isinstance(metadata, str):
        try:
            metadata = json.loads(metadata)
        except (TypeError, ValueError):
            metadata = {}
    values = set()

    def visit(value):
        if isinstance(value, dict):
            for nested in value.values():
                visit(nested)
        elif isinstance(value, (list, tuple)):
            for nested in value:
                visit(nested)
        elif value is not None:
            values.add(str(value))

    visit(metadata or {})
    return values


def _deployment_data(deployment):
    executable_path = _value(
        deployment, "executablePath", "executable_path", default=""
    )
    executable_basename = str(executable_path).rsplit("/", 1)[-1].lower()
    utility_dispatcher_entrypoint = (
        "controller" in executable_basename
        or ("epolyscat" in executable_basename and "wrapper" in executable_basename)
    )
    return {
        "id": _value(deployment, "appDeploymentId", "app_deployment_id", default=""),
        "compute_host_id": _value(
            deployment, "computeHostId", "compute_host_id", default=""
        ),
        "executable_path": executable_path,
        "controller_entrypoint": utility_dispatcher_entrypoint,
        "utility_dispatcher_entrypoint": utility_dispatcher_entrypoint,
    }


def _is_uri_collection(input_data):
    data_type = _value(input_data, "type", default=None)
    if data_type == 4:
        return True
    data_type_name = getattr(data_type, "name", "")
    return str(data_type_name or data_type).upper().endswith("URI_COLLECTION")


def audit_runtime_configuration(
    *, application_module_id, modules, interfaces, deployments, errors=None
):
    modules = list(modules or [])
    interfaces = list(interfaces or [])
    deployments = list(deployments or [])
    matching_modules = [
        module
        for module in modules
        if _value(module, "appModuleId", "app_module_id") == application_module_id
    ]
    matching_interfaces = [
        interface
        for interface in interfaces
        if application_module_id
        in (_value(interface, "applicationModules", "application_modules", default=[]) or [])
    ]
    matching_deployments = [
        deployment
        for deployment in deployments
        if _value(deployment, "appModuleId", "app_module_id") == application_module_id
    ]
    deployment_records = [
        _deployment_data(deployment) for deployment in matching_deployments
    ]
    controller_deployments = [
        deployment
        for deployment in deployment_records
        if deployment["controller_entrypoint"]
    ]

    interface = matching_interfaces[0] if len(matching_interfaces) == 1 else None
    inputs = list(
        _value(interface, "applicationInputs", "application_inputs", default=[])
        or []
    )
    inputs_by_name = {
        _value(input_data, "name", default=""): input_data for input_data in inputs
    }
    utility_selector = inputs_by_name.get("Application_Utility")
    selector_values = _metadata_values(
        _value(utility_selector, "metaData", "metadata", default={})
    )
    # Utility selectors require the dispatching controller, not the direct
    # ePolyScat executable used by module-only deployments.
    deployed = bool(controller_deployments)
    control_input_configured = any(
        input_name in inputs_by_name
        for input_name in utility_runtime_contracts.UTILITY_CONTROL_INPUT_NAMES
    )

    utilities = []
    for utility_id in EXPECTED_UTILITIES:
        contract = utility_runtime_contracts.get_utility_contract(utility_id)
        required_inputs = contract["data_input_names"]
        minimum_data_file_count = contract["minimum_data_file_count"]
        selectable = utility_id in selector_values
        direct_inputs_configured = all(
            input_name in inputs_by_name for input_name in required_inputs
        )
        generic_input_configured = (
            utility_id != "MoldenMerge"
            and any(
                input_name in inputs_by_name
                for input_name in GENERIC_POSTPROCESSING_INPUTS
            )
        )
        required_inputs_configured = (
            direct_inputs_configured or generic_input_configured
        )
        required_input_mode = (
            "direct"
            if direct_inputs_configured
            else "generic_collection"
            if generic_input_configured
            else "missing"
        )
        configured_input_names = (
            required_inputs
            if direct_inputs_configured
            else GENERIC_POSTPROCESSING_INPUTS
            if generic_input_configured
            else ()
        )
        supports_required_file_count = bool(
            required_inputs_configured
            and (
                minimum_data_file_count <= 1
                or any(
                    _is_uri_collection(inputs_by_name[input_name])
                    for input_name in configured_input_names
                    if input_name in inputs_by_name
                )
            )
        )
        ready = bool(
            matching_modules
            and interface is not None
            and deployed
            and selectable
            and control_input_configured
            and required_inputs_configured
            and supports_required_file_count
        )
        utilities.append(
            {
                "id": utility_id,
                "ready": ready,
                "selectable": selectable,
                "required_inputs": list(required_inputs),
                "required_inputs_configured": required_inputs_configured,
                "required_input_mode": required_input_mode,
                "control_input_configured": control_input_configured,
                "minimum_data_file_count": minimum_data_file_count,
                "supports_required_file_count": supports_required_file_count,
                "deployed": deployed,
            }
        )

    return {
        "ready": all(utility["ready"] for utility in utilities),
        "application_module_id": application_module_id,
        "wrapper": {
            "module_found": bool(matching_modules),
            "interface_count": len(matching_interfaces),
            "interface_id": _value(
                interface,
                "applicationInterfaceId",
                "application_interface_id",
                default="",
            ),
            "deployment_count": len(matching_deployments),
            "utility_controller_deployment_count": len(controller_deployments),
            "utility_controller_compute_host_ids": [
                deployment["compute_host_id"]
                for deployment in controller_deployments
            ],
            "deployments": deployment_records,
        },
        "utilities": utilities,
        "errors": list(errors or []),
    }
