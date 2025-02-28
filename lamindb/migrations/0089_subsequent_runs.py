from django.db import migrations, models


def update_artifact_run_relationships(apps, schema_editor):
    vendor = schema_editor.connection.vendor
    with schema_editor.connection.cursor() as cursor:
        # Step 1: Add the current run_id to the _previous_runs table if it doesn't exist
        # This preserves the latest run before we modify run_id
        cursor.execute("""
            INSERT INTO lamindb_artifact__previous_runs (artifact_id, run_id)
            SELECT a.id, a.run_id
            FROM lamindb_artifact a
            WHERE a.run_id IS NOT NULL
            AND NOT EXISTS (
                SELECT 1
                FROM lamindb_artifact__previous_runs apr
                WHERE apr.artifact_id = a.id
                AND apr.run_id = a.run_id
            );
        """)

        # Step 2: For each artifact, find the earliest run (lowest ID) and set it as the run_id
        if vendor == "sqlite":
            cursor.execute("""
                UPDATE lamindb_artifact
                SET run_id = (
                    SELECT MIN(r.id)
                    FROM lamindb_run r
                    JOIN lamindb_artifact__previous_runs apr ON r.id = apr.run_id
                    WHERE apr.artifact_id = lamindb_artifact.id
                )
                WHERE EXISTS (
                    SELECT 1
                    FROM lamindb_artifact__previous_runs apr
                    WHERE apr.artifact_id = lamindb_artifact.id
                );
            """)
        else:  # PostgreSQL
            cursor.execute("""
                UPDATE lamindb_artifact AS a
                SET run_id = subquery.min_run_id
                FROM (
                    SELECT artifact_id, MIN(run_id) as min_run_id
                    FROM lamindb_artifact__previous_runs
                    GROUP BY artifact_id
                ) AS subquery
                WHERE a.id = subquery.artifact_id
                AND EXISTS (
                    SELECT 1
                    FROM lamindb_artifact__previous_runs apr
                    WHERE apr.artifact_id = a.id
                );
            """)

        # Step 3: Remove the earliest run from the link table
        # This needs to be done after step 2 to avoid affecting the previous queries
        if vendor == "sqlite":
            cursor.execute("""
                DELETE FROM lamindb_artifact__previous_runs
                WHERE EXISTS (
                    SELECT 1 FROM lamindb_artifact a
                    WHERE lamindb_artifact__previous_runs.artifact_id = a.id
                    AND lamindb_artifact__previous_runs.run_id = a.run_id
                );
            """)
        else:  # PostgreSQL
            cursor.execute("""
                DELETE FROM lamindb_artifact__previous_runs AS apr
                USING lamindb_artifact AS a
                WHERE apr.artifact_id = a.id
                AND apr.run_id = a.run_id;
            """)


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0088_schema_components"),
    ]

    operations = [
        migrations.RunPython(
            update_artifact_run_relationships, migrations.RunPython.noop
        ),
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
                        to="lamindb.run",
                        related_name="_recreated_output_artifacts",
                        db_table="lamindb_artifact__previous_runs",  # Keep the original table name
                    ),
                ),
            ],
        ),
    ]
