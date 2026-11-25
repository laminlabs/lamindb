from django.db import migrations


def install_guard(apps, schema_editor):
    if schema_editor.connection.vendor != "postgresql":
        return
    schema_editor.execute(
        """
        CREATE FUNCTION public.lamindb_require_empty_merged_branch()
        RETURNS trigger
        LANGUAGE plpgsql
        SECURITY DEFINER
        SET search_path = pg_catalog, pg_temp
        SET row_security = off
        AS $$
        DECLARE
            source_table record;
            has_records boolean;
        BEGIN
            -- Discover membership tables, including installed schema modules.
            -- Branch annotations and provenance references are not membership.
            FOR source_table IN
                SELECT DISTINCT ns.nspname, rel.relname
                FROM pg_constraint fk
                JOIN pg_class rel ON rel.oid = fk.conrelid
                JOIN pg_namespace ns ON ns.oid = rel.relnamespace
                JOIN pg_attribute col
                    ON col.attrelid = rel.oid AND col.attnum = fk.conkey[1]
                WHERE fk.contype = 'f' AND fk.confrelid = TG_RELID
                  AND array_length(fk.conkey, 1) = 1
                  AND col.attname = 'branch_id'
                  AND EXISTS (
                      SELECT 1 FROM pg_attribute provenance
                      WHERE provenance.attrelid = rel.oid
                        AND provenance.attname = 'created_on_id'
                        AND NOT provenance.attisdropped
                  )
            LOOP
                EXECUTE format(
                    'SELECT EXISTS(SELECT 1 FROM %I.%I WHERE branch_id = $1)',
                    source_table.nspname, source_table.relname
                ) INTO has_records USING NEW.id;
                IF has_records THEN
                    RAISE EXCEPTION 'Cannot mark branch as merged while source records remain.'
                        USING ERRCODE = '23514';
                END IF;
            END LOOP;
            RETURN NEW;
        END;
        $$;
        REVOKE ALL ON FUNCTION public.lamindb_require_empty_merged_branch() FROM PUBLIC;
        CREATE TRIGGER lamindb_require_empty_merged_branch
        BEFORE UPDATE OF _status_code ON public.lamindb_branch
        FOR EACH ROW WHEN (NEW._status_code = -1)
        EXECUTE FUNCTION public.lamindb_require_empty_merged_branch();
        """,
        params=None,
    )


def uninstall_guard(apps, schema_editor):
    if schema_editor.connection.vendor == "postgresql":
        schema_editor.execute(
            "DROP TRIGGER lamindb_require_empty_merged_branch ON public.lamindb_branch;"
            "DROP FUNCTION public.lamindb_require_empty_merged_branch();"
        )


class Migration(migrations.Migration):
    dependencies = [("lamindb", "0193_v2_9_part_2")]
    operations = [migrations.RunPython(install_guard, uninstall_guard)]
