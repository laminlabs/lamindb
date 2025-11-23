import lamindb as ln
from django.db import connection


def test_simple_rename_no_parent():
    """Test renaming a record type with no parent (root level)."""

    # Create a record type
    cell_type = ln.Record(name="CellType", is_type=True).save()

    # Create a feature that uses this record type
    feature = ln.Feature(name="cell_annotation", dtype=cell_type).save()

    print(f"\nInitial feature.dtype: {repr(feature.dtype)}")

    # Check database before rename
    with connection.cursor() as cursor:
        cursor.execute(
            "SELECT id, name, dtype FROM lamindb_feature WHERE name = 'cell_annotation'"
        )
        row = cursor.fetchone()
        print(f"DB before rename: id={row[0]}, name={row[1]}, dtype={repr(row[2])}")

        # Test the conditions that will be checked in the trigger
        old_path = "CellType"  # What we expect old_path to be

        print(f"\nTesting WHERE conditions with old_path='{old_path}':")

        # Condition 1: dtype LIKE '%cat[Record[%'
        cursor.execute(
            """
            SELECT dtype LIKE %s as cond1
            FROM lamindb_feature
            WHERE name = 'cell_annotation'
        """,
            ["%cat[Record[%"],
        )
        print(f"  Condition 1 (dtype LIKE '%cat[Record[%'): {cursor.fetchone()[0]}")

        # Condition 2: dtype LIKE '%' || old_path || ']%'
        pattern2 = f"%{old_path}]%"
        cursor.execute(
            """
            SELECT dtype LIKE %s as cond2
            FROM lamindb_feature
            WHERE name = 'cell_annotation'
        """,
            [pattern2],
        )
        print(f"  Condition 2 (dtype LIKE '%{old_path}]%'): {cursor.fetchone()[0]}")

        # Combined conditions
        cursor.execute(
            """
            SELECT dtype LIKE %s AND dtype LIKE %s as both_conds
            FROM lamindb_feature
            WHERE name = 'cell_annotation'
        """,
            ["%cat[Record[%", pattern2],
        )
        print(f"  Both conditions (AND): {cursor.fetchone()[0]}")

        # Test what REPLACE would produce
        replace_pattern = f"{old_path}]"
        new_replace = "CellTypeRenamed]"
        cursor.execute(
            """
            SELECT
                dtype,
                REPLACE(dtype, %s, %s) as new_dtype
            FROM lamindb_feature
            WHERE name = 'cell_annotation'
        """,
            [replace_pattern, new_replace],
        )
        row = cursor.fetchone()
        print("\nREPLACE test:")
        print(f"  Original: {repr(row[0])}")
        print(
            f"  After REPLACE(dtype, '{replace_pattern}', '{new_replace}'): {repr(row[1])}"
        )

    # Add this right before the rename
    with connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT id, name, is_type
            FROM lamindb_record
            WHERE id = %s
        """,
            [cell_type.id],
        )
        row = cursor.fetchone()
        print(f"\nBefore rename - Record: id={row[0]}, name={row[1]}, is_type={row[2]}")

    # Rename
    cell_type.name = "CellTypeRenamed"
    cell_type.save()

    # Check after
    with connection.cursor() as cursor:
        cursor.execute(
            """
            SELECT id, name, is_type
            FROM lamindb_record
            WHERE id = %s
        """,
            [cell_type.id],
        )
        row = cursor.fetchone()
        print(f"After rename - Record: id={row[0]}, name={row[1]}, is_type={row[2]}")

        # Check the trigger's actual WHEN condition
        cursor.execute("""
            SELECT tgname, pg_get_triggerdef(oid)
            FROM pg_trigger
            WHERE tgname LIKE '%record%name%'
        """)
        for trigger in cursor.fetchall():
            print(f"\nTrigger definition:\n{trigger[1]}")

    # Check database after rename
    with connection.cursor() as cursor:
        cursor.execute(
            "SELECT id, name, dtype FROM lamindb_feature WHERE name = 'cell_annotation'"
        )
        row = cursor.fetchone()
        print(f"DB after rename: id={row[0]}, name={row[1]}, dtype={repr(row[2])}")

    # Refresh from database
    feature.refresh_from_db()

    print(f"Final feature.dtype: {repr(feature.dtype)}\n")

    with connection.cursor() as cursor:
        cursor.execute("""
            SELECT pg_get_functiondef(oid)
            FROM pg_proc
            WHERE proname = 'pgtrigger_update_feature_dtype_on_record_type_name_change_28a83'
        """)
        func_def = cursor.fetchone()
        if func_def:
            print(f"\nActual trigger function in database:\n{func_def[0]}")
        else:
            print("\nTrigger function not found!")

    assert feature.dtype == "cat[Record[CellTypeRenamed]]"
