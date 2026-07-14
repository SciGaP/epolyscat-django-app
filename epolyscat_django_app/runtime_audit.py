"""Audit the shared Airavata wrapper used by ePolyScat modules and utilities."""

import json


EXPECTED_UTILITIES = (
    "CnvMath",
    "CnvMatLab",
    "CnvLinFull",
    "MoldenMerge",
    "NRFPAD",
    "Cube2igor",
)

UTILITY_REQUIRED_INPUTS = {
    "CnvMath": ("BendOrient_Output",),
    "CnvMatLab": ("BendOrient_Output",),
    "CnvLinFull": ("DumpOut",),
    "MoldenMerge": ("molden.dat",),
    "NRFPAD": ("Cross_Section_Input_File",),
    "Cube2igor": ("Cube_Output",),
}

GENERIC_POSTPROCESSING_INPUTS = (
    "ePolyscat_Input_Data",
    "ePolyScat_Input_Data",
)


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
    return {
        "id": _value(deployment, "appDeploymentId", "app_deployment_id", default=""),
        "compute_host_id": _value(
            deployment, "computeHostId", "compute_host_id", default=""
        ),
        "executable_path": executable_path,
        "controller_entrypoint": "controller" in executable_path.lower(),
    }


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

    utilities = []
    for utility_id in EXPECTED_UTILITIES:
        required_inputs = UTILITY_REQUIRED_INPUTS[utility_id]
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
        ready = bool(
            matching_modules
            and interface is not None
            and deployed
            and selectable
            and required_inputs_configured
        )
        utilities.append(
            {
                "id": utility_id,
                "ready": ready,
                "selectable": selectable,
                "required_inputs": list(required_inputs),
                "required_inputs_configured": required_inputs_configured,
                "required_input_mode": required_input_mode,
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
