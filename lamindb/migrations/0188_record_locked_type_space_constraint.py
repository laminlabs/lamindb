from django.db import migrations

CREATE_LOCKED_TYPE_SPACE_TRIGGER = """
CREATE OR REPLACE FUNCTION enforce_record_locked_type_space() RETURNS TRIGGER AS $$
BEGIN
    IF NEW.type_id IS NOT NULL AND EXISTS (
        SELECT 1
        FROM lamindb_record r
        WHERE r.id = NEW.type_id
          AND r.is_locked
          AND r.space_id IS DISTINCT FROM NEW.space_id
    ) THEN
        RAISE EXCEPTION 'Cannot set type: record space must match locked type space';
    END IF;
    RETURN NEW;
END;
$$ LANGUAGE plpgsql;

DROP TRIGGER IF EXISTS record_locked_type_space_check ON lamindb_record;
CREATE TRIGGER record_locked_type_space_check
BEFORE INSERT OR UPDATE OF type_id, space_id
ON lamindb_record
FOR EACH ROW
WHEN (NEW.type_id IS NOT NULL)
EXECUTE FUNCTION enforce_record_locked_type_space();
"""

DROP_LOCKED_TYPE_SPACE_TRIGGER = """
DROP TRIGGER IF EXISTS record_locked_type_space_check ON lamindb_record;
DROP FUNCTION IF EXISTS enforce_record_locked_type_space();
"""


def add_locked_type_space_trigger(apps, schema_editor):
    if schema_editor.connection.vendor == "postgresql":
        schema_editor.execute(CREATE_LOCKED_TYPE_SPACE_TRIGGER)


def remove_locked_type_space_trigger(apps, schema_editor):
    if schema_editor.connection.vendor == "postgresql":
        schema_editor.execute(DROP_LOCKED_TYPE_SPACE_TRIGGER)


class Migration(migrations.Migration):
    dependencies = [
        ("lamindb", "0187_v2_4_part_2"),
    ]

    operations = [
        migrations.RunPython(
            add_locked_type_space_trigger,
            remove_locked_type_space_trigger,
        )
    ]
