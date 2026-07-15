"""Resolve the command-line contract for a registered Airavata deployment."""

import os


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


def _is_controller_entrypoint(executable_path):
    return "controller" in os.path.basename(str(executable_path or "")).lower()


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
    if _is_controller_entrypoint(executable_path):
        return {
            "exclusive": False,
            "input_names": forced_input_names,
        }

    mode = str(input_values.get("Calculation_Type") or "").upper()
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


def apply_command_line_policy(application_inputs, policy):
    input_names = set(policy.get("input_names") or ())
    exclusive = bool(policy.get("exclusive"))
    for input_data in application_inputs or ():
        input_name = _value(input_data, "name")
        if exclusive:
            input_data.requiredToAddedToCommandLine = input_name in input_names
        elif input_name in input_names:
            input_data.requiredToAddedToCommandLine = True


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
