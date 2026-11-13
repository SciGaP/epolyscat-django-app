import importlib


def _added_field_names(migration_module_name):
    migration_module = importlib.import_module(migration_module_name)
    return [
        operation.name
        for operation in migration_module.Migration.operations
        if operation.__class__.__name__ == "AddField"
    ]


def test_workflow_metadata_migration_matches_applied_schema_split():
    field_names = _added_field_names(
        "epolyscat_django_app.migrations.0024_run_workflow_metadata"
    )

    assert field_names == [
        "run_mode",
        "workflow_stage",
        "workflow_application",
        "workflow_metadata",
    ]


def test_module_and_utility_application_fields_have_followup_migration():
    migration_module = importlib.import_module(
        "epolyscat_django_app.migrations.0025_run_application_selectors"
    )
    field_names = _added_field_names(
        "epolyscat_django_app.migrations.0025_run_application_selectors"
    )

    assert migration_module.Migration.dependencies == [
        ("epolyscat_django_app", "0024_run_workflow_metadata")
    ]
    assert field_names == ["module_application", "utility_application"]


def test_workflow_child_step_fields_have_followup_migration():
    migration_module = importlib.import_module(
        "epolyscat_django_app.migrations.0026_run_workflow_child_steps"
    )
    field_names = _added_field_names(
        "epolyscat_django_app.migrations.0026_run_workflow_child_steps"
    )

    assert migration_module.Migration.dependencies == [
        ("epolyscat_django_app", "0025_run_application_selectors")
    ]
    assert field_names == [
        "parent_run",
        "workflow_step_order",
        "workflow_step_status",
    ]


def test_workflow_continuation_source_has_followup_migration():
    migration_module = importlib.import_module(
        "epolyscat_django_app.migrations.0027_run_workflow_source_run"
    )
    field_names = _added_field_names(
        "epolyscat_django_app.migrations.0027_run_workflow_source_run"
    )

    assert migration_module.Migration.dependencies == [
        ("epolyscat_django_app", "0026_run_workflow_child_steps")
    ]
    assert field_names == ["workflow_source_run"]
