# ruff: noqa: S608
from django.db import migrations, models


def update_model_run_relationships(apps, schema_editor, model_name):
    vendor = schema_editor.connection.vendor

    # Define table names based on model_name
    model_table = f"lamindb_{model_name}"
    link_table = f"lamindb_{model_name}__previous_runs"
    model_id_field = f"{model_name}_id"

    with schema_editor.connection.cursor() as cursor:
        # Step 1: Add the current run_id to the _previous_runs table if it doesn't exist
        cursor.execute(f"""
            INSERT INTO {link_table} ({model_id_field}, run_id)
            SELECT a.id, a.run_id
            FROM {model_table} a
            WHERE a.run_id IS NOT NULL
            AND NOT EXISTS (
                SELECT 1
                FROM {link_table} apr
                WHERE apr.{model_id_field} = a.id
                AND apr.run_id = a.run_id
            );
        """)

        # Step 2: For each model, find the earliest run (lowest ID) and set it as the run_id
        if vendor == "sqlite":
            cursor.execute(f"""
                UPDATE {model_table}
                SET run_id = (
                    SELECT MIN(r.id)
                    FROM lamindb_run r
                    JOIN {link_table} apr ON r.id = apr.run_id
                    WHERE apr.{model_id_field} = {model_table}.id
                )
                WHERE EXISTS (
                    SELECT 1
                    FROM {link_table} apr
                    WHERE apr.{model_id_field} = {model_table}.id
                );
            """)
        else:  # PostgreSQL
            cursor.execute(f"""
                UPDATE {model_table} AS a
                SET run_id = subquery.min_run_id
                FROM (
                    SELECT {model_id_field}, MIN(run_id) as min_run_id
                    FROM {link_table}
                    GROUP BY {model_id_field}
                ) AS subquery
                WHERE a.id = subquery.{model_id_field}
                AND EXISTS (
                    SELECT 1
                    FROM {link_table} apr
                    WHERE apr.{model_id_field} = a.id
                );
            """)

        # Step 3: Remove the earliest run from the link table
        if vendor == "sqlite":
            cursor.execute(f"""
                DELETE FROM {link_table}
                WHERE EXISTS (
                    SELECT 1 FROM {model_table} a
                    WHERE {link_table}.{model_id_field} = a.id
                    AND {link_table}.run_id = a.run_id
                );
            """)
        else:  # PostgreSQL
            cursor.execute(f"""
                DELETE FROM {link_table} AS apr
                USING {model_table} AS a
                WHERE apr.{model_id_field} = a.id
                AND apr.run_id = a.run_id;
            """)


def update_artifact_run_relationships(apps, schema_editor):
    """Migration function for artifacts."""
    update_model_run_relationships(apps, schema_editor, "artifact")


def update_collection_run_relationships(apps, schema_editor):
    """Migration function for collections."""
    update_model_run_relationships(apps, schema_editor, "collection")


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0088_schema_components"),
    ]

    operations = [
        # unrelated to subsequent runs, but related to lamindb 1.2
        # update the otype field in the artifact table
        migrations.RunSQL(
            sql="""
            UPDATE lamindb_artifact
            SET otype = 'SpatialData'
            WHERE otype = 'spatialdata';
            """
        ),
        # Migrate artifact relationships
        migrations.RunPython(
            update_artifact_run_relationships, migrations.RunPython.noop
        ),
        # Update artifact model state
        migrations.SeparateDatabaseAndState(
            # Database operations (none, to keep tables intact)
            [],
            # State operations (to update Django model only)
            [
                # Remove the old field from model state
                migrations.RemoveField(
                    model_name="artifact",
                    name="_previous_runs",
                ),
                # Add the new field with the same underlying table
                migrations.AddField(
                    model_name="artifact",
                    name="_subsequent_runs",
                    field=models.ManyToManyField(
                        "lamindb.run",
                        related_name="_recreated_artifacts",
                        db_table="lamindb_artifact__previous_runs",  # Keep the original table name
                    ),
                ),
            ],
        ),
        # Migrate collection relationships
        migrations.RunPython(
            update_collection_run_relationships, migrations.RunPython.noop
        ),
        # Update collection model state
        migrations.SeparateDatabaseAndState(
            # Database operations (none, to keep tables intact)
            [],
            # State operations (to update Django model only)
            [
                # Remove the old field from model state
                migrations.RemoveField(
                    model_name="collection",
                    name="_previous_runs",
                ),
                # Add the new field with the same underlying table
                migrations.AddField(
                    model_name="collection",
                    name="_subsequent_runs",
                    field=models.ManyToManyField(
                        "lamindb.run",
                        related_name="_recreated_collections",
                        db_table="lamindb_collection__previous_runs",  # Keep the original table name
                    ),
                ),
            ],
        ),
    ]
