from epolyscat_django_app import scientific_output_contracts


def _reader(contents_by_name):
    return lambda entry: contents_by_name.get(entry["name"])


def test_manifest_distinguishes_scheduler_artifacts_from_scientific_outputs():
    manifest = scientific_output_contracts.build_output_manifest(
        [
            {
                "name": "ePolyScat.stdout",
                "fileType": "STDOUT",
                "data-product-uri": "airavata-dp://stdout",
            },
            {
                "name": "test03dumpidy.dat",
                "fileType": "DumpOut",
                "data-product-uri": "airavata-dp://dump",
            },
        ]
    )

    assert manifest[0]["scheduler_artifact"] is True
    assert manifest[0]["scientific"] is False
    assert manifest[0]["data_product_uri"] == "airavata-dp://stdout"
    assert manifest[1]["roles"] == ["dump_out"]
    assert manifest[1]["scheduler_artifact"] is False
    assert manifest[1]["scientific"] is True


def test_gaussian_requires_normal_termination_marker():
    manifest = scientific_output_contracts.build_output_manifest(
        [{"name": "gaussian.log", "data-product-uri": "airavata-dp://log"}]
    )

    verified = scientific_output_contracts.verify_scientific_outputs(
        "Gaussian16",
        manifest,
        read_text=_reader(
            {"gaussian.log": "Normal termination of Gaussian 16 at Fri Jul 17"}
        ),
    )
    failed = scientific_output_contracts.verify_scientific_outputs(
        "Gaussian16",
        manifest,
        read_text=_reader(
            {"gaussian.log": "Error termination via Lnk1e in Gaussian"}
        ),
    )

    assert verified["status"] == "verified"
    assert failed["status"] == "failed"


def test_openmolcas_requires_success_marker_and_molden_output():
    manifest = scientific_output_contracts.build_output_manifest(
        [
            {"name": "test12.molcas.out", "data-product-uri": "airavata-dp://log"},
            {
                "name": "test12.molcas.scf.molden",
                "data-product-uri": "airavata-dp://molden",
            },
        ]
    )

    result = scientific_output_contracts.verify_scientific_outputs(
        "OpenMolcas",
        manifest,
        read_text=_reader(
            {
                "test12.molcas.out": (
                    "--- Stop Module: scf /rc=_RC_ALL_IS_WELL_ ---\nHappy landing"
                )
            }
        ),
    )

    assert result["status"] == "verified"
    assert result["evidence_files"] == [
        "test12.molcas.out",
        "test12.molcas.scf.molden",
    ]


def test_openmolcas_without_molden_is_incomplete_even_when_log_succeeds():
    manifest = scientific_output_contracts.build_output_manifest(
        [{"name": "molcas.out", "data-product-uri": "airavata-dp://log"}]
    )

    result = scientific_output_contracts.verify_scientific_outputs(
        "OpenMolcas",
        manifest,
        read_text=_reader({"molcas.out": "Happy landing"}),
    )

    assert result["status"] == "incomplete"
    assert result["reason"] == "missing_scientific_output"


def test_epolyscat_requires_completion_markers_and_scientific_result():
    manifest = scientific_output_contracts.build_output_manifest(
        [
            {
                "name": "ePolyScat.stdout",
                "fileType": "STDOUT",
                "data-product-uri": "airavata-dp://stdout",
            },
            {
                "name": "test03dumpidy.dat",
                "fileType": "DumpOut",
                "data-product-uri": "airavata-dp://dump",
            },
        ]
    )

    verified = scientific_output_contracts.verify_scientific_outputs(
        "ePolyScat",
        manifest,
        read_text=_reader({"ePolyScat.stdout": "End EDCS\nFinalize\n"}),
    )
    failed = scientific_output_contracts.verify_scientific_outputs(
        "ePolyScat",
        manifest,
        read_text=_reader({"ePolyScat.stdout": "End of input reached\n"}),
    )

    assert verified["status"] == "verified"
    assert failed["status"] == "failed"


def test_utility_requires_its_declared_output_role():
    manifest = scientific_output_contracts.build_output_manifest(
        [
            {
                "name": "converted-bendorient.dat",
                "fileType": "BendOrient",
                "data-product-uri": "airavata-dp://bend",
            }
        ]
    )

    success = scientific_output_contracts.verify_scientific_outputs(
        "CnvMath", manifest
    )
    wrong_utility = scientific_output_contracts.verify_scientific_outputs(
        "CnvLinFull", manifest
    )

    assert success["status"] == "verified"
    assert wrong_utility["status"] == "incomplete"


def test_existing_output_without_readable_log_is_unverified_not_failed():
    manifest = scientific_output_contracts.build_output_manifest(
        [{"name": "gaussian.log", "data-product-uri": "airavata-dp://log"}]
    )

    result = scientific_output_contracts.verify_scientific_outputs(
        "Gaussian16", manifest, read_text=lambda _entry: None
    )

    assert result["status"] == "unverified"
    assert result["reason"] == "evidence_unavailable"
