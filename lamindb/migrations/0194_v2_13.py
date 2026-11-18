from django.db import connection, migrations


class Migration(migrations.Migration):
    atomic = False  # CONCURRENTLY cannot run inside a transaction
    dependencies = [("lamindb", "0193_v2_9_part_2")]
    operations = []  # type: ignore


if connection.vendor == "postgresql":
    Migration.operations = [
        migrations.RunSQL(
            [
                # clears an INVALID index left by a failed build
                "DROP INDEX CONCURRENTLY IF EXISTS lamindb_artifact_tree_idx",
                """
                CREATE INDEX CONCURRENTLY lamindb_artifact_tree_idx ON lamindb_artifact
                  (space_id, storage_id, is_latest, key COLLATE "C", id)
                  INCLUDE (branch_id, uid, kind)
                """,
            ],
            reverse_sql="DROP INDEX CONCURRENTLY IF EXISTS lamindb_artifact_tree_idx",
        )
    ]
