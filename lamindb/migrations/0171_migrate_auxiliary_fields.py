"""Migrate auxiliary fields to standard Django fields."""

from django.db import migrations, models

from lamindb.models.schema import migrate_auxiliary_fields_postgres


def truncate_schema_uids(apps, schema_editor):
    """Truncate Schema UIDs from 20 to 16 chars before altering the field."""
    db_vendor = schema_editor.connection.vendor
    if db_vendor == "postgresql":
        # Raw SQL for PostgreSQL
        schema_editor.execute(
            "UPDATE lamindb_schema SET uid = LEFT(uid, 16) WHERE LENGTH(uid) > 16"
        )
    else:
        # Python fallback for SQLite
        Schema = apps.get_model("lamindb", "Schema")
        for obj in Schema.objects.filter().iterator():
            if len(obj.uid) > 16:
                obj.uid = obj.uid[:16]
                obj.save(update_fields=["uid"])


def migrate_auxiliary_fields(apps, schema_editor):
    """Migrate data from _aux['af'] to new fields."""
    db_vendor = schema_editor.connection.vendor

    if db_vendor == "postgresql":
        # Raw SQL for PostgreSQL - much faster for large datasets
        migrate_auxiliary_fields_postgres(schema_editor)
    else:
        # Python fallback for SQLite
        _migrate_auxiliary_fields_python(apps)


def _migrate_auxiliary_fields_python(apps):
    """Python-based migration for SQLite."""
    # Artifact: migrate _save_completed
    Artifact = apps.get_model("lamindb", "Artifact")
    for obj in Artifact.objects.iterator():
        af = (obj._aux or {}).get("af", {})
        obj._save_completed = af.get("0")
        if obj._aux and "af" in obj._aux:
            obj._aux.pop("af", None)
            if not obj._aux:
                obj._aux = None
        obj.save(update_fields=["_save_completed", "_aux"])

    # Run: migrate cli_args
    Run = apps.get_model("lamindb", "Run")
    for obj in Run.objects.filter(_aux__has_key="af").iterator():
        af = obj._aux.get("af", {})
        obj.cli_args = af.get("0")
        obj._aux.pop("af", None)
        if not obj._aux:
            obj._aux = None
        obj.save(update_fields=["cli_args", "_aux"])

    # Feature: migrate default_value, nullable, coerce
    Feature = apps.get_model("lamindb", "Feature")
    for obj in Feature.objects.iterator():
        af = (obj._aux or {}).get("af", {})
        if obj.is_type:
            obj.default_value = None
            obj.nullable = None
            obj.coerce = None
        else:
            obj.default_value = af.get("0")
            obj.nullable = af.get("1", True)
            obj.coerce = af.get("2", False)
        if obj._aux and "af" in obj._aux:
            obj._aux.pop("af", None)
            if not obj._aux:
                obj._aux = None
        obj.save(update_fields=["default_value", "nullable", "coerce", "_aux"])

    # Schema: migrate coerce, flexible, n_members
    Schema = apps.get_model("lamindb", "Schema")
    for obj in Schema.objects.iterator():
        af = (obj._aux or {}).get("af", {})
        if obj.is_type:
            obj.coerce = None
            obj.flexible = None
            obj.n_members = None
        else:
            obj.coerce = af.get("0")
            is_flexible_schema = obj.n_members is None or obj.n_members < 0
            if obj.n_members is not None and obj.n_members < 0:
                obj.n_members = None
            obj.flexible = af.get("2", is_flexible_schema)
        # Keep '1' (optionals) and '3' (index_feature_uid) in _aux
        if obj._aux and "af" in obj._aux:
            new_af = {}
            if "1" in af:
                new_af["1"] = af["1"]
            if "3" in af:
                new_af["3"] = af["3"]
            if new_af:
                obj._aux["af"] = new_af
            else:
                obj._aux.pop("af", None)
            if not obj._aux:
                obj._aux = None
        obj.save(update_fields=["coerce", "flexible", "n_members", "_aux"])


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0170_remove_transform_config_remove_transform_flow_and_more"),
    ]

    operations = [
        # Artifact: add _save_completed field (nullable: None=no write needed, False=in progress, True=completed)
        migrations.AddField(
            model_name="artifact",
            name="_save_completed",
            field=models.BooleanField(null=True, default=None),
        ),
        # Run: add cli_args field
        migrations.AddField(
            model_name="run",
            name="cli_args",
            field=models.CharField(max_length=1024, null=True, default=None),
        ),
        # Feature: add new fields
        migrations.AddField(
            model_name="feature",
            name="default_value",
            field=models.JSONField(null=True, default=None),
        ),
        migrations.AddField(
            model_name="feature",
            name="nullable",
            field=models.BooleanField(null=True, default=None),
        ),
        migrations.AddField(
            model_name="feature",
            name="coerce",
            field=models.BooleanField(null=True, default=None),
        ),
        # Schema: rename n to n_members and make nullable
        migrations.RenameField(
            model_name="schema",
            old_name="n",
            new_name="n_members",
        ),
        migrations.AlterField(
            model_name="schema",
            name="n_members",
            field=models.IntegerField(null=True, default=None),
        ),
        # Schema: add coerce and flexible fields
        migrations.AddField(
            model_name="schema",
            name="coerce",
            field=models.BooleanField(null=True, default=None),
        ),
        migrations.AddField(
            model_name="schema",
            name="flexible",
            field=models.BooleanField(null=True, default=None),
        ),
        # Schema: truncate UIDs before reducing max_length
        migrations.RunPython(truncate_schema_uids, migrations.RunPython.noop),
        # Schema: reduce uid max_length from 20 to 16
        migrations.AlterField(
            model_name="schema",
            name="uid",
            field=models.CharField(
                max_length=16, unique=True, db_index=True, editable=False
            ),
        ),
        # Run data migration for auxiliary fields
        migrations.RunPython(migrate_auxiliary_fields, migrations.RunPython.noop),
    ]
