from types import SimpleNamespace

from epolyscat_django_app import workflow_output_contracts


def test_scheduler_files_are_not_scientific_outputs():
    files = [
        {"name": "job_1683512072.slurm", "fileType": "STDOUT"},
        {"name": "ePolyScat.stderr", "fileType": "STDERR"},
        {"name": "ePolyScat.stdout", "fileType": "STDOUT"},
    ]

    assert all(
        workflow_output_contracts.classify_output_file(file_data) == []
        for file_data in files
    )


def test_gaussian_output_binds_to_epolyscat_convert_input():
    result = workflow_output_contracts.resolve_workflow_output_binding(
        output_files=[
            {"name": "job_1683512072.slurm"},
            {"name": "gaussian.log", "data-product-uri": "airavata-dp://log"},
            {"name": "molden.dat", "data-product-uri": "airavata-dp://molden"},
        ],
        source_application="Gaussian16",
        target_stage="ePolyScat_Run",
    )

    assert result["status"] == "ready"
    assert result["input_file_name"] == "ePolyScat_Input_Data"
    assert result["selected"]["name"] == "gaussian.log"
    assert result["data_entry_values"] == {
        "convertSource": "$pt/gaussian.log",
        "convertFormat": "gaussian",
    }


def test_openmolcas_output_binds_molden_file_to_epolyscat():
    result = workflow_output_contracts.resolve_workflow_output_binding(
        output_files=[
            {"name": "molcas.out"},
            {"name": "water.scf.molden", "data-product-uri": "airavata-dp://molden"},
        ],
        source_application="OpenMolcas",
        target_stage="ePolyScat_Run",
    )

    assert result["status"] == "ready"
    assert result["selected"]["name"] == "water.scf.molden"
    assert result["data_entry_values"]["convertSource"] == "$pt/water.scf.molden"
    assert result["data_entry_values"]["convertFormat"] == "molden"


def test_analysis_utilities_bind_only_their_scientific_input_roles():
    files = [
        {"name": "cross-sections.dat"},
        {"name": "test15dumpidy.dat", "data-product-uri": "airavata-dp://dump"},
        {"name": "bendorient.dat", "data-product-uri": "airavata-dp://bend"},
        {"name": "orientncro.out", "data-product-uri": "airavata-dp://orient"},
        {"name": "density.cube", "data-product-uri": "airavata-dp://cube"},
    ]
    expectations = {
        "CnvLinFull": ("DumpOut", "test15dumpidy.dat"),
        "CnvMath": ("BendOrient_Output", "bendorient.dat"),
        "CnvMatLab": ("BendOrient_Output", "bendorient.dat"),
        "NRFPAD": ("Cross_Section_Input_File", "orientncro.out"),
        "Cube2igor": ("Cube_Output", "density.cube"),
    }

    for application, (input_name, filename) in expectations.items():
        result = workflow_output_contracts.resolve_workflow_output_binding(
            output_files=files,
            source_application="ePolyScat",
            target_stage="Analysis",
            target_application=application,
            required_file_name=input_name,
        )
        assert result["status"] == "ready"
        assert result["input_file_name"] == input_name
        assert result["selected"]["name"] == filename


def test_missing_scientific_output_does_not_bind_generic_data_file():
    result = workflow_output_contracts.resolve_workflow_output_binding(
        output_files=[
            {"name": "cross-sections.dat"},
            {"name": "ePolyScat.stdout", "fileType": "STDOUT"},
        ],
        source_application="ePolyScat",
        target_stage="Analysis",
        target_application="CnvLinFull",
        required_file_name="DumpOut",
    )

    assert result["status"] == "missing"
    assert result["selected"] is None
    assert result["expected_role"] == "dump_out"


def test_output_classifier_accepts_object_descriptors_and_airavata_file_types():
    output = SimpleNamespace(name="result.dat", fileType="DumpOut")

    assert workflow_output_contracts.classify_output_file(output) == ["dump_out"]
