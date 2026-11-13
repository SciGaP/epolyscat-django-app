from django.db import migrations, models


def assign_tutorial_keys(apps, _schema_editor):
    """Leave legacy rows and their histories untouched.

    Schema-only on new installs. The packaged catalog uses separate v2 keys;
    databases that already ran the earlier data migration also remain intact.
    """
    pass


def clear_tutorial_keys(apps, _schema_editor):
    # A no-op reverse must not alter legacy identities on rollback.
    pass


class Migration(migrations.Migration):

    dependencies = [
        ("epolyscat_django_app", "0027_run_workflow_source_run"),
    ]

    operations = [
        migrations.AddField(
            model_name="view",
            name="tutorial_key",
            field=models.CharField(
                blank=True,
                editable=False,
                max_length=80,
                null=True,
                unique=True,
            ),
        ),
        migrations.AddField(
            model_name="run",
            name="tutorial_key",
            field=models.CharField(
                blank=True,
                editable=False,
                max_length=80,
                null=True,
                unique=True,
            ),
        ),
        migrations.RunPython(assign_tutorial_keys, clear_tutorial_keys),
    ]
