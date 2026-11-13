import json
from types import SimpleNamespace as Record
from unittest import mock

import pytest

from epolyscat_django_app.application_catalog import build_catalog, load_registry, DiscoveryUnavailable


def registry(name="ePolyScat", module_id="new-registry-id"):
    return dict(
        modules=[Record(appModuleName=name, appModuleId=module_id)],
        interfaces=[Record(applicationInterfaceId="interface", applicationModules=[module_id], applicationInputs=[])],
        deployments=[Record(appModuleId=module_id, appDeploymentId="deployment", computeHostId="cluster", executablePath="/bin/ePolyScat")],
    )


def test_ids_come_from_registry_without_gateway_specific_defaults():
    catalog = build_catalog(**registry())
    assert catalog["EPOLYSCAT"]["EPOLYSCAT_APPLICATION_ID"] == "new-registry-id"
    assert catalog["applications"] == [dict(id="ePolyScat", label="ePolyScat", kind="module", moduleId="new-registry-id", interfaceId="interface")]
    assert catalog["EPOLYSCAT"]["GAUSSIAN16_APPLICATION_ID"] is None


def test_unsupported_and_undeployed_modules_are_not_offered():
    assert build_catalog(**registry(name="Unrelated"))["applications"] == []
    data = registry()
    data["deployments"] = []
    assert build_catalog(**data)["applications"] == []


def test_ambiguous_names_require_explicit_configuration():
    data = registry()
    data["modules"].append(Record(appModuleName="ePolyScat", appModuleId="another-version"))
    assert build_catalog(**data)["applications"] == []
    catalog = build_catalog(**data, configured_ids={"EPOLYSCAT_APPLICATION_ID": "new-registry-id"})
    assert len(catalog["applications"]) == 1


def test_stale_or_inaccessible_configured_id_never_falls_back():
    catalog = build_catalog(**registry(), configured_ids={"EPOLYSCAT_APPLICATION_ID": "removed-id"})
    assert catalog["applications"] == []
    assert catalog["issues"]


def test_ambiguous_interfaces_are_not_offered():
    data = registry()
    data["interfaces"].append(Record(applicationInterfaceId="second", applicationModules=["new-registry-id"]))
    assert build_catalog(**data)["applications"] == []


@pytest.mark.parametrize("nested_config", [False, True])
def test_utility_choices_require_advertised_options_and_dispatcher(nested_config):
    data = registry()
    options = {"options": [{"value": "CnvMath"}, {"value": "Unrelated"}]}
    data["interfaces"][0].applicationInputs = [
        Record(name="Application_Utility", metaData=json.dumps({"editor": {"config": options} if nested_config else options})),
        Record(name="ePolyscat_Input_File"),
        Record(name="BendOrient_Output"),
    ]
    assert [app["id"] for app in build_catalog(**data)["applications"]] == ["ePolyScat"]
    data["deployments"][0].executablePath = "/opt/epolyscat/controller.sh"
    assert [app["id"] for app in build_catalog(**data)["applications"]] == ["ePolyScat", "CnvMath"]
    data["interfaces"][0].applicationInputs[0].metaData = '{"description": "CnvMath"}'
    assert [app["id"] for app in build_catalog(**data)["applications"]] == ["ePolyScat"]


def test_discovery_is_user_scoped_and_only_cached_on_request():
    client = Record(getAccessibleAppModules=mock.Mock(return_value=[]), getAllApplicationInterfaces=mock.Mock(return_value=[]), getAllApplicationDeployments=mock.Mock(return_value=[]))
    request = Record(airavata_client=client, authz_token="user-token")
    first = load_registry(request, "gateway")
    assert load_registry(request, "gateway") is first
    client.getAccessibleAppModules.assert_called_once_with("user-token", "gateway")
    load_registry(Record(airavata_client=client, authz_token="another-user"), "gateway")
    assert client.getAccessibleAppModules.call_count == 2


def test_api_failure_is_not_an_empty_success_or_fixed_id_fallback():
    request = Record(airavata_client=Record(getAccessibleAppModules=mock.Mock(side_effect=RuntimeError("offline"))), authz_token="token")
    with pytest.raises(DiscoveryUnavailable):
        load_registry(request, "gateway")


def test_molden_merge_is_offered_with_a_single_file_input():
    data = registry()
    data["deployments"][0].executablePath = "/opt/epolyscat/controller.sh"
    data["interfaces"][0].applicationInputs = [
        Record(name="Application_Utility", metaData={"editor": {"config": {"options": [{"value": "MoldenMerge"}]}}}),
        Record(name="ePolyscat_Input_File", type=3),
        Record(name="molden.dat", type=3),
    ]
    catalog = build_catalog(**data)
    assert "MoldenMerge" in [app["id"] for app in catalog["applications"]]
    assert not catalog["issues"]
