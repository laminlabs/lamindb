"""Tests for PostgreSQL data migrations."""

import os

import lamindb as ln
import pytest


@pytest.mark.skipif(
    os.getenv("LAMINDB_TEST_DB_VENDOR") != "postgresql",
    reason="PostgreSQL-specific migration test",
)
def test_migrate_auxiliary_fields_postgres():
    """Test PostgreSQL migration of auxiliary fields for models.

    This test verifies that migrate_auxiliary_fields_postgres correctly migrates:

    **Artifact:**
    - _save_completed from _aux['af']['0']

    **Run:**
    - cli_args from _aux['af']['0']

    **Feature:**
    - default_value from _aux['af']['0']
    - nullable from _aux['af']['1'] (default: True)
    - coerce from _aux['af']['2'] (default: False)
    - For type features, all values are set to NULL

    **Schema:**
    - coerce from _aux['af']['0']
    - flexible from _aux['af']['2'] (or computes from n_members)
    - Converts negative n_members to NULL
    - For type schemas, all values are set to NULL
    - Preserves '1' (optionals) and '3' (index_feature_uid) in _aux
    """
    from django.db import connection
    from lamindb.models.schema import migrate_auxiliary_fields_postgres

    # === Setup test data ===

    # Create a Transform and Run for testing
    transform = ln.Transform(key="test_migration_transform").save()
    run = ln.Run(transform=transform).save()

    # Create an Artifact for testing
    artifact = ln.Artifact(".gitignore", key="test_migration_artifact").save()

    # Create Features for testing (type and regular)
    type_feature = ln.Feature(
        name="TestMigrationTypeFeat", dtype=str, is_type=True
    ).save()
    regular_feature = ln.Feature(name="test_migration_regular_feat", dtype=str).save()

    # Create Schemas for testing (type and regular)
    type_schema = ln.Schema(name="TestMigrationTypeSchema", is_type=True).save()
    feature_for_schema1 = ln.Feature(
        name="test_migration_schema_feat1", dtype=str
    ).save()
    feature_for_schema2 = ln.Feature(
        name="test_migration_schema_feat2", dtype=str
    ).save()
    regular_schema = ln.Schema(
        name="TestMigrationRegularSchema",
        features=[feature_for_schema1, feature_for_schema2],
        coerce=True,
        flexible=True,
    ).save()

    # === Add _save_completed column temporarily (removed in migration 0173) ===
    with connection.cursor() as cursor:
        cursor.execute(
            """
            DO $$
            BEGIN
                IF NOT EXISTS (
                    SELECT 1 FROM information_schema.columns
                    WHERE table_name = 'lamindb_artifact' AND column_name = '_save_completed'
                ) THEN
                    ALTER TABLE lamindb_artifact ADD COLUMN _save_completed BOOLEAN;
                END IF;
            END $$;
            """
        )

    # === Set old-style _aux data to simulate pre-migration state ===
    with connection.cursor() as cursor:
        # Artifact: set _aux with af containing _save_completed value
        cursor.execute(
            """
            UPDATE lamindb_artifact
            SET _aux = '{"af": {"0": true}}'::jsonb,
                _save_completed = NULL
            WHERE id = %s
            """,
            [artifact.id],
        )

        # Run: set _aux with af containing cli_args value
        cursor.execute(
            """
            UPDATE lamindb_run
            SET _aux = '{"af": {"0": "--verbose --debug"}}'::jsonb,
                cli_args = NULL
            WHERE id = %s
            """,
            [run.id],
        )

        # Feature (type): set _aux with af keys that should result in NULL values
        cursor.execute(
            """
            UPDATE lamindb_feature
            SET _aux = '{"af": {"0": "default_val", "1": false, "2": true}}'::jsonb,
                default_value = NULL,
                nullable = NULL,
                coerce = NULL
            WHERE id = %s
            """,
            [type_feature.id],
        )

        # Feature (regular): set _aux with af keys for migration
        cursor.execute(
            """
            UPDATE lamindb_feature
            SET _aux = '{"af": {"0": "my_default", "1": false, "2": true}}'::jsonb,
                default_value = NULL,
                nullable = NULL,
                coerce = NULL
            WHERE id = %s
            """,
            [regular_feature.id],
        )

        # Schema (type): set _aux with af keys that should be cleaned
        cursor.execute(
            """
            UPDATE lamindb_schema
            SET _aux = '{"af": {"0": true, "2": false}}'::jsonb,
                coerce = NULL,
                flexible = NULL
            WHERE id = %s
            """,
            [type_schema.id],
        )

        # Schema (regular): set _aux with af keys including optionals (key "1")
        cursor.execute(
            """
            UPDATE lamindb_schema
            SET _aux = '{"af": {"0": true, "1": ["uid1", "uid2"], "2": true}}'::jsonb,
                coerce = NULL,
                flexible = NULL
            WHERE id = %s
            """,
            [regular_schema.id],
        )

    # === Run the migration function ===
    with connection.schema_editor() as schema_editor:
        migrate_auxiliary_fields_postgres(schema_editor)

    # === Refresh all objects from database ===
    run.refresh_from_db()
    type_feature.refresh_from_db()
    regular_feature.refresh_from_db()
    type_schema.refresh_from_db()
    regular_schema.refresh_from_db()

    # === Verify Artifact migration ===
    with connection.cursor() as cursor:
        cursor.execute(
            "SELECT _save_completed, _aux FROM lamindb_artifact WHERE id = %s",
            [artifact.id],
        )
        row = cursor.fetchone()
        assert row[0] is True  # _save_completed from _aux['af']['0']
        # _aux should have 'af' removed (was only key)
        assert row[1] is None or "af" not in (
            row[1] if isinstance(row[1], dict) else {}
        )

    # === Verify Run migration ===
    assert run.cli_args == "--verbose --debug"  # from _aux['af']['0']
    # _aux should have 'af' removed
    assert run._aux is None or "af" not in run._aux

    # === Verify Feature (type) migration ===
    # Type features should have all values set to NULL
    assert type_feature.default_value is None
    assert type_feature.nullable is None
    assert type_feature.coerce is None
    # _aux should have 'af' removed
    assert type_feature._aux is None or "af" not in type_feature._aux

    # === Verify Feature (regular) migration ===
    assert regular_feature.default_value == "my_default"  # from _aux['af']['0']
    assert regular_feature.nullable is False  # from _aux['af']['1']
    assert regular_feature.coerce is True  # from _aux['af']['2']
    # _aux should have 'af' removed
    assert regular_feature._aux is None or "af" not in regular_feature._aux

    # === Verify Schema (type) migration ===
    assert type_schema.coerce is None
    assert type_schema.flexible is None
    assert type_schema.n_members is None
    # _aux should either be None or not have '0' and '2' keys in 'af'
    if type_schema._aux is not None and "af" in type_schema._aux:
        assert "0" not in type_schema._aux["af"]
        assert "2" not in type_schema._aux["af"]

    # === Verify Schema (regular) migration ===
    assert regular_schema.coerce is True  # from _aux['af']['0']
    assert regular_schema.flexible is True  # from _aux['af']['2']
    # _aux should preserve key '1' (optionals)
    assert regular_schema._aux is not None
    assert "af" in regular_schema._aux
    assert "1" in regular_schema._aux["af"]
    assert regular_schema._aux["af"]["1"] == ["uid1", "uid2"]
    # Keys '0' and '2' should be removed
    assert "0" not in regular_schema._aux["af"]
    assert "2" not in regular_schema._aux["af"]

    # === Clean up: remove temporary column and delete records ===
    with connection.cursor() as cursor:
        cursor.execute(
            """
            DO $$
            BEGIN
                IF EXISTS (
                    SELECT 1 FROM information_schema.columns
                    WHERE table_name = 'lamindb_artifact' AND column_name = '_save_completed'
                ) THEN
                    ALTER TABLE lamindb_artifact DROP COLUMN _save_completed;
                END IF;
            END $$;
            """
        )

    regular_schema.delete(permanent=True)
    type_schema.delete(permanent=True)
    feature_for_schema1.delete(permanent=True)
    feature_for_schema2.delete(permanent=True)
    regular_feature.delete(permanent=True)
    type_feature.delete(permanent=True)
    artifact.delete(permanent=True)
    run.delete(permanent=True)
    transform.delete(permanent=True)
