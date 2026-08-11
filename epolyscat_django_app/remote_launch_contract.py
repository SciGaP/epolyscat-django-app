"""Resolve the command-line contract for a registered Airavata deployment."""

import os

from epolyscat_django_app import utility_runtime_contracts


EPOLYSCAT_COMMAND_INPUT_NAMES = (
    "ePolyscat_Input_File",
    "ePolyScat_Input_File",
)


def _value(instance, *names):
    for name in names:
        if isinstance(instance, dict):
            value = instance.get(name)
        else:
            value = getattr(instance, name, None)
        if value not in (None, ""):
            return value
    return ""


def _set_value(instance, value, *names):
    for name in names:
        if isinstance(instance, dict):
            if name in instance:
                instance[name] = value
                return
        elif hasattr(instance, name):
            setattr(instance, name, value)
            return
    if isinstance(instance, dict):
        instance[names[0]] = value
    else:
        setattr(instance, names[0], value)


def _is_controller_entrypoint(executable_path):
    return "controller" in os.path.basename(str(executable_path or "")).lower()


def _is_portal_wrapper_entrypoint(executable_path):
    basename = os.path.basename(str(executable_path or "")).lower()
    return basename == "epolyscat" or (
        "epolyscat" in basename and "wrapper" in basename
    )


def _utility_run_details(input_values):
    mode = str(input_values.get("Calculation_Type") or "").upper()
    if mode == "UTILITY":
        return mode, str(input_values.get("Application_Utility") or "")
    if (
        mode == "WORKFLOW"
        and input_values.get("Application_Workflow") == "Analysis"
    ):
        return mode, str(input_values.get("Application_Utility") or "")
    return mode, ""


def _is_direct_epolyscat_run(input_values):
    mode = str(input_values.get("Calculation_Type") or "").upper()
    if mode == "MODULE":
        application = str(
            input_values.get("EPOLYSCAT_Application_Module") or ""
        ).lower()
        return application == "epolyscat"
    if mode == "WORKFLOW":
        return input_values.get("Application_Workflow") == "ePolyScat_Run"
    return False


def resolve_command_line_policy(
    *, executable_path, input_values, forced_input_names=None
):
    forced_input_names = set(forced_input_names or ())
    mode, utility_id = _utility_run_details(input_values)
    if _is_controller_entrypoint(executable_path):
        if utility_id:
            return _utility_wrapper_policy(
                mode=mode,
                utility_id=utility_id,
                input_values=input_values,
            )
        return {
            "exclusive": False,
            "input_names": forced_input_names,
        }

    if utility_id and _is_portal_wrapper_entrypoint(executable_path):
        return _utility_wrapper_policy(
            mode=mode,
            utility_id=utility_id,
            input_values=input_values,
        )

    if mode == "UTILITY":
        raise ValueError(
            "A direct ePolyScat deployment does not support utility dispatch"
        )
    if not _is_direct_epolyscat_run(input_values):
        raise ValueError(
            "A direct ePolyScat deployment only supports an ePolyScat module or workflow stage"
        )

    command_inputs = {
        input_name
        for input_name in EPOLYSCAT_COMMAND_INPUT_NAMES
        if input_values.get(input_name)
    }
    if len(command_inputs) != 1:
        raise ValueError(
            "A direct ePolyScat deployment requires exactly one ePolyScat command input file"
        )
    return {
        "exclusive": True,
        "input_names": command_inputs,
    }


def _utility_wrapper_policy(*, mode, utility_id, input_values):
    utility_runtime_contracts.validate_utility_input_values(
        utility_id,
        input_values,
    )
    utility_input_order = (
        utility_runtime_contracts.resolve_utility_argument_input_order(
            utility_id,
            input_values,
        )
    )
    selector_order = ["Calculation_Type"]
    if mode == "WORKFLOW":
        selector_order.append("Application_Workflow")
    selector_order.append("Application_Utility")
    ordered_input_names = tuple(selector_order) + utility_input_order
    return {
        "exclusive": True,
        "input_names": set(ordered_input_names),
        "ordered_input_names": ordered_input_names,
    }


def apply_command_line_policy(application_inputs, policy):
    input_names = set(policy.get("input_names") or ())
    exclusive = bool(policy.get("exclusive"))
    input_order = {
        input_name: position
        for position, input_name in enumerate(
            policy.get("ordered_input_names") or ()
        )
    }
    for input_data in application_inputs or ():
        input_name = _value(input_data, "name")
        if exclusive:
            _set_value(
                input_data,
                input_name in input_names,
                "required_to_added_to_command_line",
                "requiredToAddedToCommandLine",
            )
        elif input_name in input_names:
            _set_value(
                input_data,
                True,
                "required_to_added_to_command_line",
                "requiredToAddedToCommandLine",
            )
        if input_name in input_order:
            _set_value(
                input_data,
                input_order[input_name],
                "input_order",
                "inputOrder",
            )


def find_deployment_executable_path(
    deployments, *, application_module_id, compute_resource_id
):
    matches = [
        deployment
        for deployment in deployments or ()
        if _value(deployment, "appModuleId", "app_module_id")
        == application_module_id
        and _value(deployment, "computeHostId", "compute_host_id")
        == compute_resource_id
    ]
    if not matches:
        raise ValueError(
            "No application deployment is registered for the selected module and compute resource"
        )
    if len(matches) > 1:
        raise ValueError(
            "Multiple application deployments are registered for the selected module and compute resource"
        )

    executable_path = _value(
        matches[0], "executablePath", "executable_path"
    )
    if not executable_path:
        raise ValueError("The selected application deployment has no executable path")
    return str(executable_path)
