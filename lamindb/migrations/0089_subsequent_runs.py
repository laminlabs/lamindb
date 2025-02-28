from django.db import migrations, models


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0088_schema_components"),
    ]

    operations = [
        migrations.RunSQL(
            """
            -- For each artifact, find the earliest run (lowest ID) and set it as the run_id
            UPDATE lamindb_artifact a
            SET run_id = (
                SELECT MIN(r.id)
                FROM lamindb_run r
                JOIN lamindb_artifact__previous_runs apr ON r.id = apr.run_id
                WHERE apr.artifact_id = a.id
            )
            WHERE EXISTS (
                SELECT 1
                FROM lamindb_artifact__previous_runs apr
                WHERE apr.artifact_id = a.id
            );

            -- Add the current run_id to the _previous_runs table if it doesn't exist
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
            """,
            """
            -- No rollback needed for setting run_id as it's based on existing data
            -- If needed, we could delete the newly added run relationships:
            DELETE FROM lamindb_artifact__previous_runs
            WHERE (artifact_id, run_id) IN (
                SELECT a.id, a.run_id
                FROM lamindb_artifact a
                WHERE a.run_id IS NOT NULL
            );
            """,
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
