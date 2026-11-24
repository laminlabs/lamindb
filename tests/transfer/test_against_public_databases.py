"""Transfer from public LaminDB instances. These do not need authentication."""

from datetime import UTC, datetime
from pathlib import Path

import lamindb as ln
from lamindb.models.feature import parse_dtype


def _describe_source(connection):
    details = (
        f"Schema transfer source alias={connection.alias} backend={connection.vendor} "
        f"lamin_user={ln.setup.settings.user.handle}"
    )
    try:
        if connection.vendor == "sqlite":
            path = Path(connection.settings_dict["NAME"])
            modified = datetime.fromtimestamp(path.stat().st_mtime, UTC)
            return f"{details} file={path.name} modified={modified.isoformat()}"
        if connection.vendor == "postgresql":
            with connection.cursor() as cursor:
                cursor.execute("SELECT current_database(), current_user")
                database, role = cursor.fetchone()
            return f"{details} database={database} role={role}"
    except Exception as exc:
        return f"{details} metadata_error={type(exc).__name__}"
    return details


def test_schema_transfer_defaults_to_annotations(connected_bionty, capsys):
    # Switching instances in the fixtures reloads Django's connection registry.
    from django.core.exceptions import ObjectDoesNotExist
    from django.db import connections

    db = ln.DB("laminlabs/lamindata")
    schema_uid = "pnQvQVcQ417bfmVq"
    schemas = db.Schema
    source = connections[schemas.db]
    description = _describe_source(source)
    with capsys.disabled():
        print(description, flush=True)
    try:
        remote_schema = schemas.get(schema_uid)
    except ObjectDoesNotExist as exc:
        try:
            with source.cursor() as cursor:
                cursor.execute(
                    "SELECT EXISTS (SELECT 1 FROM lamindb_schema WHERE uid = %s)",
                    [schema_uid],
                )
                exists = bool(cursor.fetchone()[0])
            exc.add_note(
                f"Direct lookup on {source.alias} for {schema_uid}: exists={exists}"
            )
        except Exception as diagnostic_error:
            exc.add_note(f"Direct lookup failed with {type(diagnostic_error).__name__}")
        raise
    remote_member_names = remote_schema.members.to_list("name")
    assert len(remote_member_names) > 0

    remote_perturbation_dtype = remote_schema.members.get(
        name="perturbation"
    )._dtype_str
    remote_perturbation_type_uid = parse_dtype(remote_perturbation_dtype)[0]["type_uid"]
    remote_perturbation_type_name = db.Record.get(remote_perturbation_type_uid).name

    existing_local = ln.Schema.filter(uid=schema_uid).one_or_none()
    if existing_local is not None:
        existing_local.delete(permanent=True)

    existing_local_perturbation_type = ln.Record.filter(
        name=remote_perturbation_type_name, is_type=True
    ).one_or_none()
    if existing_local_perturbation_type is not None:
        ln.Record.filter(type=existing_local_perturbation_type).delete(permanent=True)
        existing_local_perturbation_type.delete(permanent=True)

    assert (
        ln.Record.filter(name=remote_perturbation_type_name, is_type=True).count() == 0
    )
    assert ln.Record.filter(type__name=remote_perturbation_type_name).count() == 0

    record_only = db.Schema.get(schema_uid).save(transfer="record")
    assert record_only.members.count() == 0

    transferred = db.Schema.get(schema_uid).save()
    transferred_member_names = transferred.members.to_list("name")
    assert transferred_member_names == remote_member_names

    perturbation_feature = transferred.members.get(name="perturbation")
    assert perturbation_feature.dtype_as_str.startswith("cat[Record[")
    assert perturbation_feature.dtype_as_object is not None

    # The dtype only needs the Perturbation type. Its records stay behind.
    assert (
        ln.Record.filter(name=remote_perturbation_type_name, is_type=True).count() == 1
    )
    assert ln.Record.filter(type__name=remote_perturbation_type_name).count() == 0

    before_count = transferred.links_feature.count()
    transferred_repeat = db.Schema.get(schema_uid).save()
    assert transferred_repeat.id == transferred.id
    assert transferred_repeat.links_feature.count() == before_count


def test_transfer_from_lamin_site_assets(connected_bionty):
    # An instance with fewer modules than the local one.
    artifact = ln.Artifact.connect("laminlabs/lamin-site-assets").get(
        "lgRNHNtMxjU0y8nIagt7"
    )
    # Implicit saving through load() also covers cache() and open() for run inputs.
    artifact.load()
    # storage=False: the file stays in the source instance.
    # Permanent deletion fails here because of the schema module mismatch.
    artifact.delete(storage=False)
