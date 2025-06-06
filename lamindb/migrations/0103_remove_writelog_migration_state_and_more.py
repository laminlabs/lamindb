# Generated by Django 5.2 on 2025-05-29 12:02

from django.db import migrations


def fix_artifact_kind(apps, schema_editor):
    Artifact = apps.get_model("lamindb", "Artifact")
    Artifact.objects.filter(kind="__lamindb__").update(kind="__lamindb_run__")


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0102_remove_writelog_branch_code_and_more"),
    ]

    operations = [
        migrations.RunPython(fix_artifact_kind),
        migrations.RemoveField(
            model_name="writelog",
            name="migration_state",
        ),
        migrations.RemoveField(
            model_name="writelog",
            name="table",
        ),
        migrations.RemoveField(
            model_name="writelog",
            name="branch",
        ),
        migrations.RemoveField(
            model_name="writelog",
            name="space",
        ),
        migrations.DeleteModel(
            name="WriteLogLock",
        ),
        migrations.DeleteModel(
            name="MigrationState",
        ),
        migrations.DeleteModel(
            name="TableState",
        ),
        migrations.DeleteModel(
            name="WriteLog",
        ),
    ]
