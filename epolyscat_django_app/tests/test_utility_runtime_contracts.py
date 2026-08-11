import pytest

from epolyscat_django_app import utility_runtime_contracts


def test_utility_contracts_match_the_manual_data_protocols():
    expected = {
        "CnvMath": ("BendOrient_Output", 1),
        "CnvMatLab": ("BendOrient_Output", 1),
        "CnvLinFull": ("DumpOut", 1),
        "MoldenMerge": ("molden.dat", 2),
        "NRFPAD": ("Cross_Section_Input_File", 1),
        "Cube2igor": ("Cube_Output", 1),
    }

    assert tuple(utility_runtime_contracts.UTILITY_CONTRACTS) == tuple(expected)
    for utility_id, (data_input_name, minimum_file_count) in expected.items():
        contract = utility_runtime_contracts.get_utility_contract(utility_id)

        assert contract["control_input_names"] == (
            "ePolyscat_Input_File",
            "ePolyScat_Input_File",
        )
        assert contract["data_input_names"][0] == data_input_name
        assert contract["minimum_data_file_count"] == minimum_file_count
        assert contract["executable_name"].lower().startswith(utility_id.lower())


def test_utility_input_validation_requires_control_and_scientific_data_files():
    with pytest.raises(ValueError, match="control file"):
        utility_runtime_contracts.validate_utility_input_values(
            "CnvMath",
            {"BendOrient_Output": "airavata-dp://bend"},
        )

    with pytest.raises(ValueError, match="BendOrient_Output"):
        utility_runtime_contracts.validate_utility_input_values(
            "CnvMath",
            {"ePolyscat_Input_File": "airavata-dp://control"},
        )


def test_molden_merge_requires_two_distinct_molden_files():
    with pytest.raises(ValueError, match="at least 2"):
        utility_runtime_contracts.validate_utility_input_values(
            "MoldenMerge",
            {
                "ePolyscat_Input_File": "airavata-dp://control",
                "molden.dat": "airavata-dp://first",
            },
        )

    utility_runtime_contracts.validate_utility_input_values(
        "MoldenMerge",
        {
            "ePolyscat_Input_File": "airavata-dp://control",
            "molden.dat": "airavata-dp://first,airavata-dp://second",
        },
    )


def test_utility_input_names_prefer_specific_data_without_duplicate_generic_aliases():
    input_values = {
        "Calculation_Type": "UTILITY",
        "Application_Utility": "CnvLinFull",
        "ePolyscat_Input_File": "airavata-dp://control",
        "DumpOut": "airavata-dp://dump",
        "ePolyscat_Input_Data": "airavata-dp://dump",
        "ePolyScat_Input_Data": "airavata-dp://dump",
    }

    assert utility_runtime_contracts.resolve_utility_argument_input_names(
        "CnvLinFull",
        input_values,
    ) == {"ePolyscat_Input_File", "DumpOut"}


def test_utility_input_names_accept_generic_collection_for_legacy_runs():
    input_values = {
        "ePolyScat_Input_File": "airavata-dp://control",
        "ePolyScat_Input_Data": "airavata-dp://first,airavata-dp://second",
    }

    assert utility_runtime_contracts.resolve_utility_argument_input_names(
        "MoldenMerge",
        input_values,
    ) == {"ePolyScat_Input_File", "ePolyScat_Input_Data"}


def test_unknown_utility_is_rejected():
    with pytest.raises(ValueError, match="Unsupported ePolyScat utility"):
        utility_runtime_contracts.get_utility_contract("UnknownUtility")
