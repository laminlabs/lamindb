"""Migrate auxiliary fields to standard Django fields."""

from django.db import migrations, models


def migrate_auxiliary_fields(apps, schema_editor):
    """Migrate data from _aux['af'] to new fields and truncate Schema UIDs."""
    # Artifact: migrate _is_saved_to_storage_location -> _save_completed
    Artifact = apps.get_model("lamindb", "Artifact")
    for obj in Artifact.objects.iterator():
        af = (obj._aux or {}).get("af", {})
        # Set from _aux or default to False
        obj._save_completed = af.get("0", False)
        # Clean up _aux
        if obj._aux and "af" in obj._aux:
            obj._aux.pop("af", None)
            if not obj._aux:
                obj._aux = None
        obj.save(update_fields=["_save_completed", "_aux"])

    # Run: migrate cli_args
    Run = apps.get_model("lamindb", "Run")
    for obj in Run.objects.filter(_aux__has_key="af").iterator():
        af = obj._aux.get("af", {})
        obj.cli_args = af.get("0")  # None if not present
        obj._aux.pop("af", None)
        if not obj._aux:
            obj._aux = None
        obj.save(update_fields=["cli_args", "_aux"])

    # Feature: migrate default_value, nullable, coerce
    Feature = apps.get_model("lamindb", "Feature")
    for obj in Feature.objects.iterator():
        af = (obj._aux or {}).get("af", {})

        if obj.is_type:
            # Type-like features: set all to NULL
            obj.default_value = None
            obj.nullable = None
            obj.coerce = None
        else:
            # Regular features: migrate values or use defaults
            obj.default_value = af.get("0")  # None if not present
            obj.nullable = af.get("1", True)  # Default True
            obj.coerce = af.get("2", False)  # Default False

        # Clean up _aux
        if obj._aux and "af" in obj._aux:
            obj._aux.pop("af", None)
            if not obj._aux:
                obj._aux = None

        obj.save(update_fields=["default_value", "nullable", "coerce", "_aux"])

    # Schema: migrate coerce, flexible, set n_members to NULL for types, truncate UID
    # Note: n has already been renamed to n_members by the schema migration
    Schema = apps.get_model("lamindb", "Schema")
    for obj in Schema.objects.iterator():
        af = (obj._aux or {}).get("af", {})

        # Truncate UID if longer than 16 chars
        if len(obj.uid) > 16:
            obj.uid = obj.uid[:16]

        if obj.is_type:
            # Type-like schemas: set all to NULL
            obj.coerce = None
            obj.flexible = None
            obj.n_members = None
        else:
            # Regular schemas: migrate values or use defaults
            obj.coerce = af.get("0", False)  # Default False
            # flexible default is n_members < 0 for backward compat
            obj.flexible = af.get("2", obj.n_members < 0 if obj.n_members else False)
            # n_members already has the value from the renamed n field

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

        obj.save(update_fields=["uid", "coerce", "flexible", "n_members", "_aux"])


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0170_remove_transform_config_remove_transform_flow_and_more"),
    ]

    operations = [
        # Artifact: add _save_completed field
        migrations.AddField(
            model_name="artifact",
            name="_save_completed",
            field=models.BooleanField(default=False),
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
        # Schema: reduce uid max_length from 20 to 16
        migrations.AlterField(
            model_name="schema",
            name="uid",
            field=models.CharField(
                max_length=16, unique=True, db_index=True, editable=False
            ),
        ),
        # Run data migration
        migrations.RunPython(migrate_auxiliary_fields, migrations.RunPython.noop),
    ]
