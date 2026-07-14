import json
from types import SimpleNamespace

from epolyscat_django_app import runtime_audit


def _input(name, options=None):
    metadata = {}
    if options is not None:
        metadata = {"editor": {"options": [{"value": value} for value in options]}}
    return SimpleNamespace(name=name, metaData=json.dumps(metadata))


def test_runtime_audit_reports_all_remote_utility_contracts_ready():
    module = SimpleNamespace(appModuleId="epolyscat-module", appModuleName="ePolyScat")
    interface = SimpleNamespace(
        applicationInterfaceId="epolyscat-interface",
        applicationModules=["epolyscat-module"],
        applicationInputs=[
            _input("Calculation_Type", ["MODULE", "UTILITY", "WORKFLOW"]),
            _input("Application_Utility", runtime_audit.EXPECTED_UTILITIES),
            _input("BendOrient_Output"),
            _input("DumpOut"),
            _input("molden.dat"),
            _input("Cross_Section_Input_File"),
            _input("Cube_Output"),
        ],
    )
    deployment = SimpleNamespace(
        appDeploymentId="deployment-id",
        appModuleId="epolyscat-module",
        computeHostId="expanse",
        executablePath="/opt/epolyscat/controller.sh",
    )

    result = runtime_audit.audit_runtime_configuration(
        application_module_id="epolyscat-module",
        modules=[module],
        interfaces=[interface],
        deployments=[deployment],
    )

    assert result["ready"] is True
    assert result["wrapper"]["module_found"] is True
    assert result["wrapper"]["interface_count"] == 1
    assert result["wrapper"]["deployment_count"] == 1
    assert result["wrapper"]["utility_controller_deployment_count"] == 1
    assert result["wrapper"]["deployments"][0]["controller_entrypoint"] is True
    assert [utility["id"] for utility in result["utilities"]] == list(
        runtime_audit.EXPECTED_UTILITIES
    )
    assert all(utility["ready"] for utility in result["utilities"])


def test_runtime_audit_reports_missing_selector_and_required_input_separately():
    interface = SimpleNamespace(
        applicationInterfaceId="epolyscat-interface",
        applicationModules=["epolyscat-module"],
        applicationInputs=[
            _input("Application_Utility", ["CnvMath"]),
            _input("BendOrient_Output"),
        ],
    )

    result = runtime_audit.audit_runtime_configuration(
        application_module_id="epolyscat-module",
        modules=[],
        interfaces=[interface],
        deployments=[],
    )
    utilities = {utility["id"]: utility for utility in result["utilities"]}

    assert result["ready"] is False
    assert utilities["CnvMath"]["selectable"] is True
    assert utilities["CnvMath"]["required_inputs_configured"] is True
    assert utilities["CnvMath"]["deployed"] is False
    assert utilities["CnvLinFull"]["selectable"] is False
    assert utilities["CnvLinFull"]["required_inputs_configured"] is False


def test_runtime_audit_accepts_generic_postprocessing_input_collection():
    module = SimpleNamespace(appModuleId="epolyscat-module")
    interface = SimpleNamespace(
        applicationInterfaceId="epolyscat-interface",
        applicationModules=["epolyscat-module"],
        applicationInputs=[
            _input("Application_Utility", runtime_audit.EXPECTED_UTILITIES),
            _input("ePolyscat_Input_Data"),
            _input("molden.dat"),
        ],
    )
    deployment = SimpleNamespace(
        appDeploymentId="deployment-id",
        appModuleId="epolyscat-module",
        computeHostId="expanse",
        executablePath="/opt/epolyscat/controller.sh",
    )

    result = runtime_audit.audit_runtime_configuration(
        application_module_id="epolyscat-module",
        modules=[module],
        interfaces=[interface],
        deployments=[deployment],
    )
    utilities = {utility["id"]: utility for utility in result["utilities"]}

    assert result["ready"] is True
    assert utilities["CnvLinFull"]["required_input_mode"] == "generic_collection"
    assert utilities["MoldenMerge"]["required_input_mode"] == "direct"


def test_runtime_audit_does_not_treat_direct_epolyscat_binary_as_utility_dispatcher():
    module = SimpleNamespace(appModuleId="epolyscat-module")
    interface = SimpleNamespace(
        applicationInterfaceId="epolyscat-interface",
        applicationModules=["epolyscat-module"],
        applicationInputs=[
            _input("Application_Utility", runtime_audit.EXPECTED_UTILITIES),
            _input("ePolyscat_Input_Data"),
            _input("molden.dat"),
        ],
    )
    deployment = SimpleNamespace(
        appDeploymentId="deployment-id",
        appModuleId="epolyscat-module",
        computeHostId="frontera",
        executablePath="/opt/epolyscat/ePolyScat",
    )

    result = runtime_audit.audit_runtime_configuration(
        application_module_id="epolyscat-module",
        modules=[module],
        interfaces=[interface],
        deployments=[deployment],
    )

    assert result["ready"] is False
    assert result["wrapper"]["deployment_count"] == 1
    assert result["wrapper"]["utility_controller_deployment_count"] == 0
    assert all(utility["deployed"] is False for utility in result["utilities"])
