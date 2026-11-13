"""Packaged ePolyScat tutorial discovery and shared-run support."""

import json

from pathlib import Path

from django.db import transaction
from django.db.models import Prefetch

from epolyscat_django_app import models, output_presentation_contracts
from epolyscat_django_app.tutorial_contract import TUTORIAL_VIEW_KEY, TUTORIAL_RUN_KEY_PREFIX


PACKAGE_DIR = Path(__file__).resolve().parent
MANUAL_TUTORIALS_DIR = PACKAGE_DIR / "data" / "epolyscat" / "manual_tutorials"
MANIFEST_PATH = MANUAL_TUTORIALS_DIR / "manifest.json"
EXPECTED_CASE_IDS = ("test01", "test04", "test08", "test12", "test36")
ARTIFACT_FIELDS = (
    "script",
    "runtime_data",
    "standard_output",
    "reference_input",
)
REQUIRED_CATALOG_FIELDS = (
    "id",
    "name",
    "description",
    "input_format",
    *ARTIFACT_FIELDS,
    "convert_format",
)
ALLOWED_CONVERT_FORMATS = frozenset(("molden", "gaussian", "gamess"))
CATALOG_REVISION = 1


def load_catalog():
    """Load the curated manual sample catalog in display order."""

    with MANIFEST_PATH.open(encoding="utf-8") as manifest_stream:
        manifest = json.load(manifest_stream)
    if not isinstance(manifest, dict) or manifest.get("version") != 1:
        raise ValueError("Tutorial catalog must use manifest version 1")
    cases = list(manifest["cases"])
    if tuple(case.get("id") for case in cases) != EXPECTED_CASE_IDS:
        raise ValueError("Tutorial catalog must contain exactly the approved cases")

    normalized_cases = []
    catalog_root = MANUAL_TUTORIALS_DIR.resolve()
    for case in cases:
        missing_fields = [
            field for field in REQUIRED_CATALOG_FIELDS if not case.get(field)
        ]
        if missing_fields:
            raise ValueError(
                "Tutorial catalog entry is missing required fields: "
                + ", ".join(missing_fields)
            )
        if case["convert_format"] not in ALLOWED_CONVERT_FORMATS:
            raise ValueError("Tutorial convert format is not supported")

        normalized_case = dict(case)
        for field in ARTIFACT_FIELDS:
            filename = str(case.get(field) or "")
            if not filename or Path(filename).name != filename:
                raise ValueError("Tutorial artifact paths must be direct basenames")
        case_path = (catalog_root / normalized_case["id"]).resolve()
        if case_path.parent != catalog_root:
            raise ValueError("Tutorial case paths must be direct children")
        for field in ARTIFACT_FIELDS:
            artifact_path = case_path / normalized_case[field]
            if not artifact_path.is_file():
                raise ValueError(
                    f"Tutorial declared artifact is missing: {normalized_case[field]}"
                )
            if artifact_path.is_symlink() or artifact_path.resolve().parent != case_path:
                raise ValueError(
                    "Tutorial declared artifact must remain inside its case directory"
                )
        normalized_case["case_path"] = case_path
        normalized_cases.append(normalized_case)
    return normalized_cases


def discover_cases():
    """Return the validated case directories in manifest display order."""

    return [case["case_path"] for case in load_catalog()]


def _tutorial_run_key(case_id):
    return f"{TUTORIAL_RUN_KEY_PREFIX}{case_id}"


def _workflow_metadata(case):
    return {
        "tutorialStandalone": True,
        "tutorialCatalogRevision": CATALOG_REVISION,
        "tutorialCaseId": case["id"],
        "tutorialInputFormat": case["input_format"],
        "tutorialConvertFormat": case["convert_format"],
        "tutorialRuntimeDataName": case["runtime_data"],
        "tutorialReferenceFiles": [case["reference_input"]],
        "output_files": [case["standard_output"]],
        "file_groups": [
            {
                "label": "Inputs",
                "files": [case["script"], case["runtime_data"]],
            },
            {"label": "Output", "files": [case["standard_output"]]},
            {"label": "Reference", "files": [case["reference_input"]]},
        ],
        "selected_file": case["script"],
    }


def _tutorial_view_queryset(for_update=False):
    queryset = models.View.objects
    if for_update:
        queryset = queryset.select_for_update()
    return queryset.prefetch_related(
        Prefetch(
            "runs",
            queryset=models.Run.objects.prefetch_related(
                Prefetch(
                    "inputs",
                    queryset=models.Input.objects.prefetch_related("files"),
                )
            ),
        )
    )


def _matching_inputs(run, input_type, name):
    return [
        input_instance
        for input_instance in run.inputs.all()
        if input_instance.type == input_type and input_instance.name == name
    ]


def _run_inputs_are_current(run, case):
    expected_input_identities = {
        ("radio input", "Calculation_Type"),
        ("radio input", "EPOLYSCAT_Application_Module"),
        ("files", "ePolyscat_Input_File"),
        ("files", "ePolyScat_Input_Data"),
    }
    inputs = list(run.inputs.all())
    if (
        len(inputs) != len(expected_input_identities)
        or {(input_instance.type, input_instance.name) for input_instance in inputs}
        != expected_input_identities
    ):
        return False

    expected_selectors = {
        "Calculation_Type": "MODULE",
        "EPOLYSCAT_Application_Module": "ePolyScat",
    }
    for name, value in expected_selectors.items():
        matches = _matching_inputs(run, "radio input", name)
        if len(matches) != 1 or matches[0].value != value:
            return False

    expected_files = {
        "ePolyscat_Input_File": case["script"],
        "ePolyScat_Input_Data": case["runtime_data"],
    }
    file_inputs = [
        input_instance
        for input_instance in run.inputs.all()
        if input_instance.type == "files"
    ]
    if (
        len(file_inputs) != len(expected_files)
        or {input_instance.name for input_instance in file_inputs}
        != set(expected_files)
    ):
        return False
    for input_name, filename in expected_files.items():
        matches = _matching_inputs(run, "files", input_name)
        if len(matches) != 1:
            return False
        files = list(matches[0].files.all())
        if len(files) != 1:
            return False
        if files[0].name != filename or files[0].data_product_uri:
            return False
    return True


def _run_is_current(run, case):
    relative_path = case["case_path"].relative_to(PACKAGE_DIR).as_posix()
    expected_fields = {
        "tutorial_key": _tutorial_run_key(case["id"]),
        "owner_id": None,
        "root_id": None,
        "experiment_id": None,
        "parent_run_id": None,
        "workflow_source_run_id": None,
        "airavata_project_id": None,
        "inpc_data_product_uri": None,
        "directory": "",
        "name": case["name"],
        "number": case["id"],
        "description": case["description"],
        "filepath": relative_path,
        "deleted": False,
        "run_mode": "module",
        "module_application": "ePolyScat",
        "workflow_stage": "",
        "workflow_application": "",
        "utility_application": "",
        "workflow_metadata": _workflow_metadata(case),
    }
    return all(
        getattr(run, field) == expected
        for field, expected in expected_fields.items()
    ) and not run.executions.exists() and _run_inputs_are_current(run, case)


def _tutorial_view_is_current(tutorial_view, cases):
    if (
        tutorial_view.tutorial_key != TUTORIAL_VIEW_KEY
        or tutorial_view.owner_id is not None
        or tutorial_view.type != "tutorial"
        or tutorial_view.name != "Tutorials"
        or tutorial_view.deleted
    ):
        return False

    runs = list(tutorial_view.runs.all())
    runs_by_key = {run.tutorial_key: run for run in runs}
    expected_keys = {_tutorial_run_key(case["id"]) for case in cases}
    if len(runs) != len(cases) or set(runs_by_key) != expected_keys:
        return False
    return all(
        _run_is_current(runs_by_key[_tutorial_run_key(case["id"])], case)
        for case in cases
    )


def _sync_selector_input(run, name, value):
    matches = list(
        models.Input.objects.filter(
            run=run,
            type="radio input",
            name=name,
        ).order_by("pk")
    )
    input_instance = matches[0] if matches else models.Input(
        run=run,
        type="radio input",
        name=name,
    )
    if input_instance.value != value or input_instance.pk is None:
        input_instance.value = value
        input_instance.save()
    if len(matches) > 1:
        models.Input.objects.filter(pk__in=[item.pk for item in matches[1:]]).delete()
    return input_instance


def _sync_file_input(run, input_name, filename):
    matches = list(
        models.Input.objects.filter(
            run=run,
            type="files",
            name=input_name,
        ).order_by("pk")
    )
    input_instance = matches[0] if matches else models.Input(
        run=run,
        type="files",
        name=input_name,
        value=None,
    )
    if input_instance.pk is None or input_instance.value is not None:
        input_instance.value = None
        input_instance.save()
    if len(matches) > 1:
        models.Input.objects.filter(pk__in=[item.pk for item in matches[1:]]).delete()
    input_instance.files.exclude(name=filename).delete()
    file_instance = input_instance.files.filter(name=filename).order_by("pk").first()
    if file_instance is None:
        file_instance = models.File(input=input_instance, name=filename)
    file_instance.data_product_uri = ""
    file_instance.save()
    input_instance.files.filter(name=filename).exclude(pk=file_instance.pk).delete()
    return input_instance


def ensure_tutorial_view():
    """Synchronize and return the shared Tutorials view."""

    cases = load_catalog()
    tutorial_view = _tutorial_view_queryset().filter(
        tutorial_key=TUTORIAL_VIEW_KEY
    ).first()
    if tutorial_view is not None and _tutorial_view_is_current(
        tutorial_view,
        cases,
    ):
        return models.View.objects.get(pk=tutorial_view.pk)

    with transaction.atomic():
        tutorial_view = (
            models.View.objects.select_for_update()
            .filter(tutorial_key=TUTORIAL_VIEW_KEY)
            .first()
        )
        if tutorial_view is None:
            # Never adopt a legacy view: its membership/history may be personal.
            tutorial_view, _ = models.View.objects.get_or_create(
                tutorial_key=TUTORIAL_VIEW_KEY,
                defaults={"name": "Tutorials", "type": "tutorial", "owner": None},
            )

        tutorial_view = _tutorial_view_queryset(for_update=True).get(
            pk=tutorial_view.pk
        )
        if _tutorial_view_is_current(tutorial_view, cases):
            return models.View.objects.get(pk=tutorial_view.pk)

        if (
            tutorial_view.name != "Tutorials"
            or tutorial_view.type != "tutorial"
            or tutorial_view.owner_id is not None
            or tutorial_view.deleted
        ):
            tutorial_view.name = "Tutorials"
            tutorial_view.type = "tutorial"
            tutorial_view.owner = None
            tutorial_view.deleted = False
            tutorial_view.save(
                update_fields=["name", "type", "owner", "deleted", "updated"]
            )

        synchronized_runs = []
        for case in cases:
            relative_path = case["case_path"].relative_to(PACKAGE_DIR).as_posix()
            tutorial_key = _tutorial_run_key(case["id"])
            run = (
                models.Run.objects.select_for_update()
                .filter(tutorial_key=tutorial_key)
                .first()
            )
            if run is None:
                # v2 keys deliberately do not reuse v1 or unkeyed legacy records.
                run, _ = models.Run.objects.get_or_create(
                    tutorial_key=tutorial_key,
                    defaults={"owner": None, "filepath": relative_path},
                )
            if (
                run.owner_id is not None or run.root_id is not None
                or run.experiment_id is not None or run.parent_run_id is not None
                or run.workflow_source_run_id is not None
                or run.airavata_project_id is not None
                or run.inpc_data_product_uri is not None or run.directory
                or run.executions.exists()
            ):
                raise ValueError("Packaged tutorial source has personal or execution associations")

            run.owner = None
            run.filepath = relative_path
            run.name = case["name"]
            run.number = case["id"]
            run.description = case["description"]
            run.deleted = False
            run.run_mode = "module"
            run.module_application = "ePolyScat"
            run.workflow_stage = ""
            run.workflow_application = ""
            run.utility_application = ""
            run.workflow_metadata = _workflow_metadata(case)
            run.save()

            canonical_inputs = [
                _sync_selector_input(run, "Calculation_Type", "MODULE"),
                _sync_selector_input(
                    run,
                    "EPOLYSCAT_Application_Module",
                    "ePolyScat",
                ),
                _sync_file_input(run, "ePolyscat_Input_File", case["script"]),
                _sync_file_input(run, "ePolyScat_Input_Data", case["runtime_data"]),
            ]
            run.inputs.exclude(
                pk__in=[input_instance.pk for input_instance in canonical_inputs]
            ).delete()
            synchronized_runs.append(run)

        tutorial_view.runs.set(synchronized_runs)
        return models.View.objects.get(pk=tutorial_view.pk)


def case_for_run(run):
    """Resolve a synchronized tutorial run to its validated manifest entry."""

    if run.owner_id is not None or not run.views.filter(
        type="tutorial",
        owner=None,
        tutorial_key=TUTORIAL_VIEW_KEY,
    ).exists():
        raise FileNotFoundError("Run is not a packaged tutorial")

    case_id = (run.workflow_metadata or {}).get("tutorialCaseId")
    if run.tutorial_key != _tutorial_run_key(case_id):
        raise FileNotFoundError("Run is not a canonical packaged tutorial")
    try:
        catalog = load_catalog()
    except (OSError, KeyError, TypeError, ValueError) as exc:
        raise FileNotFoundError("Tutorial catalog is unavailable") from exc
    case = next(
        (catalog_case for catalog_case in catalog if catalog_case["id"] == case_id),
        None,
    )
    if case is None:
        raise FileNotFoundError("Tutorial run is not in the packaged catalog")

    candidate = (PACKAGE_DIR / run.filepath).resolve()
    if candidate != case["case_path"]:
        raise FileNotFoundError("Tutorial case is outside the packaged catalog")
    return case


def case_path_for_run(run):
    """Resolve a synchronized tutorial run to its packaged case directory."""

    return case_for_run(run)["case_path"]


def input_file_names(run):
    case = case_for_run(run)
    return [case["script"], case["runtime_data"]]


def output_file_names(run):
    return [case_for_run(run)["standard_output"]]


def reference_file_names(run):
    return [case_for_run(run)["reference_input"]]


def _file_descriptor(filename):
    return {
        "name": filename,
        "filename": filename,
        "data_product_uri": "",
        "data-product-uri": "",
        "viewable": output_presentation_contracts.is_viewable_file(filename),
    }


def input_file_descriptors(run):
    return [_file_descriptor(filename) for filename in input_file_names(run)]


def output_file_descriptors(run):
    return output_presentation_contracts.annotate_output_files(
        [_file_descriptor(filename) for filename in output_file_names(run)]
    )


def reference_file_descriptors(run):
    return [_file_descriptor(filename) for filename in reference_file_names(run)]


def open_case_file(run, filename):
    """Open one packaged tutorial file after strict basename validation."""

    requested_name = str(filename or "")
    if not requested_name or Path(requested_name).name != requested_name:
        raise FileNotFoundError("Invalid tutorial filename")

    case = case_for_run(run)
    case_path = case["case_path"]
    allowed_names = {case[field] for field in ARTIFACT_FIELDS}
    if requested_name not in allowed_names:
        raise FileNotFoundError("Tutorial file is not declared in the catalog")
    file_path = (case_path / requested_name).resolve()
    if file_path.parent != case_path or not file_path.is_file():
        raise FileNotFoundError("Tutorial file does not exist")
    return file_path.open("rb")
