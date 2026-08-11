from django.db import migrations, models
import django.db.models.deletion


class Migration(migrations.Migration):

    dependencies = [
        ("epolyscat_django_app", "0026_run_workflow_child_steps"),
    ]

    operations = [
        migrations.AddField(
            model_name="run",
            name="workflow_source_run",
            field=models.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.PROTECT,
                related_name="workflow_continuations",
                to="epolyscat_django_app.run",
            ),
        ),
    ]
