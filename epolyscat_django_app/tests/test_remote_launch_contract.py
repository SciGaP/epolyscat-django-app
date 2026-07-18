from types import SimpleNamespace

import pytest

from epolyscat_django_app import remote_launch_contract


def test_controller_deployment_keeps_interface_command_line_contract():
    policy = remote_launch_contract.resolve_command_line_policy(
        executable_path="/opt/epolyscat/bin/epolyscat_controller.sh",
        input_values={
            "Calculation_Type": "MODULE",
            "EPOLYSCAT_Application_Module": "ePolyScat",
            "ePolyscat_Input_File": "airavata-dp://commands",
        },
        forced_input_names={"Gaussian_Inputs"},
    )

    assert policy == {
        "exclusive": False,
        "input_names": {"Gaussian_Inputs"},
    }


@pytest.mark.parametrize(
    "input_values",
    [
        {
            "Calculation_Type": "MODULE",
            "EPOLYSCAT_Application_Module": "ePolyScat",
            "ePolyScat_Input_Data": "airavata-dp://molden",
            "ePolyscat_Input_File": "airavata-dp://commands",
        },
        {
            "Calculation_Type": "WORKFLOW",
            "Application_Workflow": "ePolyScat_Run",
            "ePolyScat_Input_Data": "airavata-dp://molden",
            "ePolyscat_Input_File": "airavata-dp://commands",
        },
    ],
)
def test_direct_epolyscat_deployment_passes_only_command_input(input_values):
    policy = remote_launch_contract.resolve_command_line_policy(
        executable_path="/opt/epolyscat/bin/ePolyScat",
        input_values=input_values,
        forced_input_names={"Calculation_Type", "ePolyScat_Input_Data"},
    )

    assert policy == {
        "exclusive": True,
        "input_names": {"ePolyscat_Input_File"},
    }


def test_exclusive_policy_disables_selectors_and_data_file_arguments():
    application_inputs = [
        SimpleNamespace(
            name="Calculation_Type", requiredToAddedToCommandLine=True
        ),
        SimpleNamespace(
            name="EPOLYSCAT_Application_Module",
            requiredToAddedToCommandLine=True,
        ),
        SimpleNamespace(
            name="ePolyScat_Input_Data", requiredToAddedToCommandLine=True
        ),
        SimpleNamespace(
            name="ePolyscat_Input_File", requiredToAddedToCommandLine=True
        ),
    ]

    remote_launch_contract.apply_command_line_policy(
        application_inputs,
        {
            "exclusive": True,
            "input_names": {"ePolyscat_Input_File"},
        },
    )

    assert {
        input_data.name: input_data.requiredToAddedToCommandLine
        for input_data in application_inputs
    } == {
        "Calculation_Type": False,
        "EPOLYSCAT_Application_Module": False,
        "ePolyScat_Input_Data": False,
        "ePolyscat_Input_File": True,
    }


def test_nonexclusive_policy_preserves_interface_flags_and_forces_named_inputs():
    application_inputs = [
        SimpleNamespace(
            name="Calculation_Type", requiredToAddedToCommandLine=True
        ),
        SimpleNamespace(
            name="Gaussian_Inputs", requiredToAddedToCommandLine=False
        ),
    ]

    remote_launch_contract.apply_command_line_policy(
        application_inputs,
        {
            "exclusive": False,
            "input_names": {"Gaussian_Inputs"},
        },
    )

    assert application_inputs[0].requiredToAddedToCommandLine is True
    assert application_inputs[1].requiredToAddedToCommandLine is True


def test_direct_epolyscat_deployment_rejects_utility_dispatch():
    with pytest.raises(ValueError, match="does not support utility dispatch"):
        remote_launch_contract.resolve_command_line_policy(
            executable_path="/opt/epolyscat/bin/ePolyScat",
            input_values={
                "Calculation_Type": "UTILITY",
                "Application_Utility": "CnvMath",
            },
        )


def test_frontera_wrapper_passes_only_the_selected_utility_contract():
    policy = remote_launch_contract.resolve_command_line_policy(
        executable_path="/opt/epolyscat/bin/epolyscat_wrapper.sh",
        input_values={
            "Calculation_Type": "UTILITY",
            "Application_Utility": "CnvMath",
            "ePolyscat_Input_File": "airavata-dp://control",
            "BendOrient_Output": "airavata-dp://bend",
            "ePolyscat_Input_Data": "airavata-dp://bend",
        },
    )

    assert policy == {
        "exclusive": True,
        "input_names": {
            "Calculation_Type",
            "Application_Utility",
            "ePolyscat_Input_File",
            "BendOrient_Output",
        },
    }


def test_frontera_wrapper_supports_workflow_analysis_utility_contract():
    policy = remote_launch_contract.resolve_command_line_policy(
        executable_path="/opt/epolyscat/bin/epolyscat_wrapper.sh",
        input_values={
            "Calculation_Type": "WORKFLOW",
            "Application_Workflow": "Analysis",
            "Application_Utility": "CnvLinFull",
            "ePolyscat_Input_File": "airavata-dp://control",
            "DumpOut": "airavata-dp://dump",
        },
    )

    assert policy == {
        "exclusive": True,
        "input_names": {
            "Calculation_Type",
            "Application_Workflow",
            "Application_Utility",
            "ePolyscat_Input_File",
            "DumpOut",
        },
    }


def test_frontera_wrapper_rejects_incomplete_utility_contract():
    with pytest.raises(ValueError, match="control file"):
        remote_launch_contract.resolve_command_line_policy(
            executable_path="/opt/epolyscat/bin/epolyscat_wrapper.sh",
            input_values={
                "Calculation_Type": "UTILITY",
                "Application_Utility": "Cube2igor",
                "Cube_Output": "airavata-dp://cube",
            },
        )


def test_deployment_lookup_uses_module_and_compute_resource():
    deployments = [
        SimpleNamespace(
            appModuleId="epolyscat-module",
            computeHostId="expanse",
            executablePath="/opt/epolyscat/bin/epolyscat_controller.sh",
        ),
        SimpleNamespace(
            appModuleId="epolyscat-module",
            computeHostId="frontera",
            executablePath="/opt/epolyscat/bin/ePolyScat",
        ),
        SimpleNamespace(
            appModuleId="other-module",
            computeHostId="frontera",
            executablePath="/opt/other/bin/run",
        ),
    ]

    assert remote_launch_contract.find_deployment_executable_path(
        deployments,
        application_module_id="epolyscat-module",
        compute_resource_id="frontera",
    ) == "/opt/epolyscat/bin/ePolyScat"


def test_deployment_lookup_rejects_missing_or_ambiguous_registration():
    deployment = SimpleNamespace(
        appModuleId="epolyscat-module",
        computeHostId="frontera",
        executablePath="/opt/epolyscat/bin/ePolyScat",
    )

    with pytest.raises(ValueError, match="No application deployment"):
        remote_launch_contract.find_deployment_executable_path(
            [],
            application_module_id="epolyscat-module",
            compute_resource_id="frontera",
        )

    with pytest.raises(ValueError, match="Multiple application deployments"):
        remote_launch_contract.find_deployment_executable_path(
            [deployment, deployment],
            application_module_id="epolyscat-module",
            compute_resource_id="frontera",
        )
