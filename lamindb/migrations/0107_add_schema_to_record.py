# Generated by Django 5.2 on 2025-06-30 12:42

import django.db.models.deletion
from django.db import migrations

import lamindb.base.fields


def migrate_sheets_to_records(apps, schema_editor):
    """Migrate Sheet records to Record table as type records."""
    with schema_editor.connection.cursor() as cursor:
        # Insert sheets as records with is_type=True
        cursor.execute("""
            INSERT INTO lamindb_record (uid, name, description, schema_id, is_type, created_by_id, created_at, updated_at, run_id)
            SELECT uid, name, description, schema_id, TRUE, created_by_id, created_at, updated_at, run_id
            FROM lamindb_sheet;
        """)

        # Update records that were linked to sheets to link to the new record types
        cursor.execute("""
            UPDATE lamindb_record
            SET type_id = (
                SELECT r.id
                FROM lamindb_record r
                JOIN lamindb_sheet s ON r.uid = s.uid
                WHERE s.id = lamindb_record.sheet_id
            )
            WHERE sheet_id IS NOT NULL;
        """)


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0106_transfer_data_migration"),
    ]

    operations = [
        migrations.RemoveConstraint(
            model_name="record",
            name="unique_name",
        ),
        migrations.AddField(
            model_name="record",
            name="schema",
            field=lamindb.base.fields.ForeignKey(
                blank=True,
                null=True,
                on_delete=django.db.models.deletion.CASCADE,
                related_name="records",
                to="lamindb.schema",
            ),
        ),
        migrations.AlterField(
            model_name="record",
            name="is_type",
            field=lamindb.base.fields.BooleanField(
                blank=True, db_index=True, default=False
            ),
        ),
        migrations.RunPython(
            migrate_sheets_to_records,
        ),
        migrations.AddField(
            model_name="record",
            name="_sort_order",
            field=django.db.models.FloatField(null=True, default=None),
        ),
    ]
