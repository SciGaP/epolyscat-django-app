from epolyscat_django_app import output_presentation_contracts


def test_file_name_declarations_preserve_arbitrary_output_names_and_last_binding():
    declarations = output_presentation_contracts.parse_file_name_declarations(
        """
        # Output names are selected by the scientist.
        FileName 'PlotData' '$po/cross sections.dat' 'REWIND'
        FileName 'PlotData2D' '$po/angular-map.dat'
        FileName 'PlotData' '$po/rebound.dat' # later PlotData binding
        """
    )

    assert declarations == {
        "cross sections.dat": "PlotData",
        "angular-map.dat": "PlotData2D",
        "rebound.dat": "PlotData",
    }


def test_output_contract_uses_declared_file_type_instead_of_filename_guessing():
    annotated = output_presentation_contracts.annotate_output_files(
        [
            {
                "name": "cross sections.dat",
                "data-product-uri": "airavata-dp://one-dimensional",
            },
            {
                "name": "angular-map.dat",
                "data-product-uri": "airavata-dp://two-dimensional",
            },
            {
                "name": "matrix.idy",
                "fileType": "MatrixElements",
                "data-product-uri": "airavata-dp://matrix",
            },
            {
                "name": "job.slurm",
                "data-product-uri": "airavata-dp://scheduler",
            },
        ],
        declared_file_types={
            "cross sections.dat": "PlotData",
            "angular-map.dat": "PlotData2D",
        },
    )

    one_dimensional, two_dimensional, matrix, scheduler = annotated
    assert one_dimensional["file_type"] == "PlotData"
    assert one_dimensional["semantic_roles"] == ["plot_data_1d"]
    assert one_dimensional["plottable"] is True
    assert one_dimensional["plot_contract"] == {
        "dimension": 1,
        "x_axis": "0",
        "y_axes": "1",
        "flags": "-linY",
    }

    assert two_dimensional["file_type"] == "PlotData2D"
    assert two_dimensional["semantic_roles"] == ["plot_data_2d"]
    assert two_dimensional["plottable"] is True
    assert two_dimensional["plot_contract"] == {
        "dimension": 2,
        "x_axis": "0,1",
        "y_axes": "2",
        "flags": "",
    }

    assert matrix["viewable"] is True
    assert matrix["plottable"] is False
    assert scheduler["viewable"] is True
    assert scheduler["plottable"] is False


def test_only_manual_outputs_supported_by_generic_plotter_are_plottable():
    annotated = output_presentation_contracts.annotate_output_files(
        [
            {"name": "mfdcs.dat", "fileType": "MFDCS"},
            {"name": "mfdcs-full.dat", "fileType": "MFDCSFull"},
            {"name": "geometry.dat", "fileType": "MFDCSGeom"},
        ]
    )

    assert annotated[0]["plottable"] is True
    assert annotated[0]["plot_contract"]["dimension"] == 2
    assert annotated[1]["plottable"] is False
    assert annotated[2]["plottable"] is False
