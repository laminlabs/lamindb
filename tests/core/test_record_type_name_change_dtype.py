import lamindb as ln
from django.db import connection


def test_simple_rename_no_parent():
    """Test renaming a record type with no parent (root level)."""

    with connection.cursor() as cursor:
        cursor.execute("SET client_min_messages TO NOTICE;")

    # Create a record type
    cell_type = ln.Record(name="CellType", is_type=True).save()

    # Create a feature that uses this record type
    feature = ln.Feature(name="cell_annotation", dtype=cell_type).save()

    # Rename the record type
    cell_type.name = "CellTypeRenamed"
    cell_type.save()

    # Refresh from database
    feature.refresh_from_db()

    # Assert the dtype was updated
    assert feature.dtype == "cat[Record[CellTypeRenamed]]"
