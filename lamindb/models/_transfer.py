from __future__ import annotations

from typing import TYPE_CHECKING, Any

import lamindb_setup as ln_setup
from django.db import ProgrammingError
from django.db.models import QuerySet as DjangoQuerySet
from lamin_utils import logger
from lamindb_setup._connect_instance import get_owner_name_from_identifier

from ..errors import NoWriteAccess, ValidationError
from .sqlrecord import BaseSQLRecord, Space, SQLRecord

if TYPE_CHECKING:
    from .run import Run

REGISTRY_UNIQUE_FIELD = {"storage": "root", "ulabel": "name"}


def update_fk_to_default_db(
    records: SQLRecord | list[SQLRecord] | DjangoQuerySet,
    fk: str,
    using: str | None,
    transfer_logs: dict,
    *,
    transfer_annotations: bool = True,
):
    # here in case it is an iterable, we are checking only a single record
    # and set the same fks for all other records because we do this only
    # for certain fks where they have to the same for the whole bulk
    # see transfer_fk_to_default_db_bulk
    # todo: but this has to be changed i think, it is not safe as it is now - Sergei
    record = records[0] if isinstance(records, (list, DjangoQuerySet)) else records
    if getattr(record, f"{fk}_id", None) is not None:
        # Map the source space by uid. Do not substitute the current space:
        # that would change who can access the object.
        if fk == "space":
            source_space = getattr(record, fk)
            fk_record_default = Space.filter(uid=source_space.uid).one_or_none()
            if fk_record_default is None:
                obj = (
                    f"{record.__class__.__name__}(uid={record.uid!r})"
                    if getattr(record, "uid", None)
                    else record.__class__.__name__
                )
                target = ln_setup.settings.instance.slug
                raise NoWriteAccess(
                    f"Could not map space {source_space.name!r} of object {obj}.\n"
                    f"Please attach space {source_space.name!r} to the target database {target!r}."
                )
        # process non-space fks
        else:
            fk_record = getattr(record, fk)
            if fk in {"created_by", "schema", "type"}:
                print(
                    f"transfer {type(record).__name__} {getattr(record, 'uid', None)} "
                    f".{fk} → {type(fk_record).__name__} {getattr(fk_record, 'uid', None)} "
                    f"{getattr(fk_record, 'handle', None) or getattr(fk_record, 'name', '')}",
                    flush=True,
                )
            field = REGISTRY_UNIQUE_FIELD.get(fk, "uid")
            pre_existing_fk_record_default = fk_record.__class__.filter(
                **{field: getattr(fk_record, field)}
            ).one_or_none()
            # A Record is only valid in a type that is already on the target,
            # because that type carries a schema. Every other missing type,
            # including a ULabel type, is a stub.
            if fk == "type" and pre_existing_fk_record_default is None:
                is_data_record = record.__class__.__name__ == "Record" and not getattr(
                    record, "is_type", False
                )
                if is_data_record:
                    type_name = getattr(fk_record, "name", None) or fk_record.uid
                    type_uid = getattr(fk_record, "uid", None)
                    raise ValueError(
                        f"Please transfer type {type_name!r} first: "
                        f"{fk_record.__class__.__name__}(uid={type_uid!r})"
                    )
                from copy import copy

                pre_existing_fk_record_default = transfer_to_default_db(
                    copy(fk_record),
                    using,
                    transfer_logs=transfer_logs,
                    stub=True,
                )
            from copy import copy

            fk_record_default = copy(fk_record)
            # A schema FK is part of the row. Its members are annotations.
            if fk_record.__class__.__name__ == "Schema" and transfer_annotations:
                from .schema import transfer_schema_with_members

                fk_record_default = transfer_schema_with_members(
                    fk_record_default, using, transfer_logs=transfer_logs
                )
            elif pre_existing_fk_record_default is None:
                transfer_to_default_db(
                    fk_record_default,
                    using,
                    save=True,
                    transfer_logs=transfer_logs,
                    transfer_annotations=transfer_annotations,
                )
            else:
                fk_record_default = pre_existing_fk_record_default
        # re-set the fks to the newly saved ones in the default db
        if isinstance(records, (list, DjangoQuerySet)):
            for r in records:
                setattr(r, f"{fk}", None)
                setattr(r, f"{fk}_id", fk_record_default.id)
        else:
            setattr(records, f"{fk}", None)
            setattr(records, f"{fk}_id", fk_record_default.id)


FKBULK = [
    "organism",
    "source",
    "report",  # Run
]


def transfer_fk_to_default_db_bulk(
    records: list | DjangoQuerySet, using: str | None, transfer_logs: dict
):
    for fk in FKBULK:
        update_fk_to_default_db(records, fk, using, transfer_logs=transfer_logs)


def get_transfer_run(record) -> Run:
    from lamindb import settings
    from lamindb.core._context import context
    from lamindb.models import Run, Transform
    from lamindb.models._lineage import WARNING_RUN_TRANSFORM

    slug = record._state.db
    owner, name = get_owner_name_from_identifier(slug)
    cache_using_filepath = (
        ln_setup.settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
    )
    if not cache_using_filepath.exists():
        raise SystemExit("Need to call .connect() before")
    instance_uid = cache_using_filepath.read_text().split("\n")[0]
    key = f"__lamindb_transfer__/{instance_uid}"
    uid = instance_uid + "0000"
    transform = Transform.filter(uid=uid).one_or_none()
    if transform is None:
        search_names = settings.creation.search_names
        settings.creation.search_names = False
        # TODO: consider renaming to "Sync from"
        transform = Transform(  # type: ignore
            uid=uid, description=f"Transfer from `{slug}`", key=key, kind="function"
        ).save()
        settings.creation.search_names = search_names
    # use the global run context to get the initiated_by_run run id
    if context.run is not None:
        initiated_by_run = context.run
    else:
        if not settings.creation.artifact_silence_missing_run_warning:
            logger.warning(WARNING_RUN_TRANSFORM)
        initiated_by_run = None
    # it doesn't seem to make sense to create new runs for every transfer
    run = Run.filter(transform=transform, initiated_by_run=initiated_by_run).first()
    if run is None:
        run = Run(
            transform=transform,
            initiated_by_run=initiated_by_run,
            status="started",
        ).save()  # type: ignore
        run.initiated_by_run = initiated_by_run  # so that it's available in memory
    return run


TRANSFER_MODES = {"sqlrecord", "notes", "annotations"}


def normalize_transfer_config(
    transfer_config: str | None, *, default_annotations: bool = False
) -> str:
    """Map transfer= to sqlrecord | notes | annotations.

    ``transfer="record"`` is kept as an alias for ``sqlrecord`` until LaminDB v3.
    Schema still defaults to ``annotations`` when ``transfer`` is omitted.
    """
    if transfer_config is None:
        return "annotations" if default_annotations else "sqlrecord"
    if transfer_config == "record":
        logger.warning(
            "transfer='record' is deprecated; use transfer='sqlrecord'. "
            "The alias will be removed in LaminDB v3."
        )
        return "sqlrecord"
    if transfer_config not in TRANSFER_MODES:
        raise ValueError(
            "transfer should be one of 'sqlrecord', 'notes', 'annotations' "
            f"(or deprecated 'record'), not {transfer_config!r}"
        )
    return transfer_config


def transfer_notes(record_on_default, source_db, source_pk) -> None:
    """Copy the latest readme block from the source SQLRecord."""
    if source_pk is None or not hasattr(record_on_default, "ablocks"):
        return
    source = record_on_default.__class__.objects.using(source_db).get(pk=source_pk)
    src_block = source.ablocks.filter(kind="readme", is_latest=True).first()
    if src_block is None or not src_block.content:
        return
    existing = record_on_default.ablocks.filter(kind="readme", is_latest=True).first()
    if existing is not None and existing.content == src_block.content:
        return
    fk_name = record_on_default.__class__.__name__.lower()
    record_on_default.ablocks.model(
        **{fk_name: record_on_default, "kind": "readme", "content": src_block.content}
    ).save()


def _user_annotation_field(feature) -> str:
    """Field `_add_values` looks up for a User feature. Defaults to handle."""
    if feature is None:
        return "handle"
    from .feature import parse_dtype

    return parse_dtype(feature._dtype_str)[0]["field_str"]


def _user_registry_write_forbidden(error: Exception) -> bool:
    message = str(error).lower()
    return "row-level security" in message or "permission denied" in message


def _save_transferred_record(record):
    """Insert a transferred row.

    `User` inserts are allowed like any other table. Until every instance
    has that policy, a rejected insert still means the person must be added
    as a collaborator so their `User` row exists, then the sync re-run.
    """
    try:
        record.save()
    except ProgrammingError as error:
        if record.__class__.__name__ != "User" or not _user_registry_write_forbidden(
            error
        ):
            raise
        handle = record.handle
        raise NoWriteAccess(
            f"Cannot write user {handle!r} (uid {record.uid!r}) to the target User registry.\n"
            f"Make {handle!r} a collaborator on this instance so they get a User entry, then re-run the sync."
        ) from None


def _map_user_annotation(source_user, feature, transfer_logs: dict):
    """Map a source User annotation onto the target User registry by uid."""
    from copy import copy

    saved = transfer_to_default_db(
        copy(source_user),
        None,
        save=True,
        transfer_logs=transfer_logs,
    )
    local = saved if saved is not None else type(source_user).get(uid=source_user.uid)
    return getattr(local, _user_annotation_field(feature))


def _linked_feature_values(record) -> list[tuple[Any, Any]]:
    """Feature values from link rows, keyed by the feature row rather than its name.

    Names are not unique. Several categorical links for one feature are one list.
    A JSON list stays one value because it is stored as a single JSON cell.
    """
    grouped: dict[int, list] = {}
    features: dict[int, Any] = {}
    for rel in record._meta.related_objects:
        accessor = rel.get_accessor_name()
        if not accessor or not str(accessor).startswith("values_"):
            continue
        for link in getattr(record, accessor).all():
            feature = link.feature
            features[feature.id] = feature
            grouped.setdefault(feature.id, []).append(link.value)
    return [
        (features[feature_id], vals[0] if len(vals) == 1 else vals)
        for feature_id, vals in grouped.items()
    ]


def transfer_record_feature_values(
    record_on_default, source_db, source_pk, using, transfer_logs
):
    from copy import copy

    from .feature import Feature, parse_dtype

    if source_pk is None:
        return
    source = record_on_default.__class__.objects.using(source_db).get(pk=source_pk)
    linked_values = _linked_feature_values(source)
    if not linked_values:
        return

    def _transfer_entity(value, feature=None):
        if type(value).__name__ == "User":
            # User is BaseSQLRecord, not SQLRecord. Return the feature field
            # (handle by default) so _add_values can look the user up.
            return _map_user_annotation(value, feature, transfer_logs)
        # A linked record is a stub (uid, name, type, created_by). Transferring
        # that record later fills its remaining fields.
        if type(value).__name__ in {"Record", "ULabel"}:
            from copy import copy

            return transfer_to_default_db(
                copy(value),
                using,
                transfer_logs=transfer_logs,
                stub=True,
            )
        return value.save(transfer="annotations")

    def _prepare(value, feature=None):
        # Link rows are reloaded from the database. Categorical values are
        # related records; JSON values are lists and scalars, not ndarrays.
        if isinstance(value, (list, tuple, set, frozenset)):
            return [
                prepared
                for v in value
                if (prepared := _prepare(v, feature)) is not None
            ]
        if isinstance(value, (SQLRecord, BaseSQLRecord)) and value._state.db not in (
            None,
            "default",
        ):
            return _transfer_entity(value, feature)
        return value

    prepared_by_uid = {}
    feature_objects = []
    for src_feature, value in linked_values:
        transferred = transfer_to_default_db(
            copy(src_feature), using, save=True, transfer_logs=transfer_logs
        )
        local_feature = (
            transferred if transferred is not None else Feature.get(uid=src_feature.uid)
        )
        dtype = local_feature._dtype_str or ""
        if dtype.startswith("cat") or dtype.startswith("list[cat"):
            try:
                parse_dtype(dtype)
            except ValidationError as err:
                raise ValueError(
                    f"cannot transfer feature {local_feature.uid!r} ({dtype}): "
                    "the target instance does not have the required schema module loaded "
                    "(e.g. run: lamin settings modules set bionty). "
                    'Pass transfer="sqlrecord" to sync the object without annotations.'
                ) from err
        prepared_by_uid[local_feature.uid] = _prepare(value, local_feature)
        feature_objects.append(local_feature)

    # Do not run ExperimentalDictCurator: source values are already valid, and
    # the target session may not have every module (e.g. bionty) imported.
    # set so _add_values can `del` it (normal set_values always sets this attr)
    record_on_default._mapped_feature_update_fields = set()
    record_on_default.features._remove_values()
    record_on_default.features._add_values(
        feature_objects,
        {},
        values_by_feature_uid=prepared_by_uid,
    )


_STUB_FKS = {"created_by", "type", "space"}


def transfer_to_default_db(
    record: SQLRecord,
    using: str | None,
    *,
    transfer_logs: dict,
    save: bool = False,
    transfer_fk: bool = True,
    transfer_annotations: bool = True,
    stub: bool = False,
) -> SQLRecord | None:
    if record._state.db is None or record._state.db == "default":
        return None
    # Dtype text is not a foreign key. Follow it even when this feature row
    # is already on the target, so a re-transfer picks up schema__uid refs.
    if record.__class__.__name__ == "Feature":
        from .feature import transfer_feature_dtypes

        transfer_feature_dtypes(record, using, transfer_logs=transfer_logs)
    registry = record.__class__
    logger.debug(f"transferring {registry.__name__} record {record.uid} to default db")
    record_on_default = registry.objects.filter(uid=record.uid).one_or_none()
    record_str = f"{record.__class__.__name__}(uid='{record.uid}')"
    if transfer_logs["run"] is None:
        transfer_logs["run"] = get_transfer_run(record)
    # A link keeps the row that is already there. Transferring the record
    # itself writes the source fields onto that row.
    filling = (
        record_on_default is not None
        and not stub
        and record.__class__.__name__ in {"Record", "ULabel"}
    )
    if record_on_default is not None and not filling:
        transfer_logs["mapped"].append(record_str)
        return record_on_default
    if filling:
        from copy import copy

        print(f"transfer fill stub {record_str}", flush=True)
        record = copy(record)
    else:
        transfer_logs["transferred"].append(record_str)
        if stub:
            print(
                f"transfer stub {record_str} {getattr(record, 'name', '')}",
                flush=True,
            )

    # run & transform stay on the transfer run; created_by is transferred
    # like any other foreign key, including the User row it points at.
    run = transfer_logs["run"]
    if hasattr(record, "run_id"):
        record.run = None
        record.run_id = run.id
    # deal with denormalized transform FK on artifact and collection
    if hasattr(record, "transform_id"):
        record.transform = None
        record.transform_id = run.transform_id
    fk_fields = [
        i.name
        for i in record._meta.fields
        if i.get_internal_type() == "ForeignKey"
        if i.name not in {"run", "transform", "branch"}
    ]
    if not transfer_fk:
        # don't transfer fk fields that are already bulk transferred
        fk_fields = [fk for fk in fk_fields if fk not in FKBULK]
    if stub:
        # Identity, type, and creator. The remaining fields are filled when
        # this record is transferred itself.
        fk_fields = [fk for fk in fk_fields if fk in _STUB_FKS]
        for name in ("description", "reference", "reference_type", "schema_id"):
            if hasattr(record, name):
                setattr(record, name, None)
    for fk in fk_fields:
        update_fk_to_default_db(
            record,
            fk,
            using,
            transfer_logs=transfer_logs,
            transfer_annotations=transfer_annotations,
        )
    # FK ids were remapped to the default DB; drop tracked *_id originals so save
    # logic does not treat remapping as a user-requested field change.
    if (original_values := getattr(record, "_original_values", None)) is not None:
        for key in [key for key in original_values if key.endswith("_id")]:
            del original_values[key]
    record._state.db = "default"
    if filling:
        for field in record._meta.concrete_fields:
            if field.primary_key:
                continue
            setattr(record_on_default, field.attname, getattr(record, field.attname))
        _save_transferred_record(record_on_default)
        return record_on_default
    record.id = None
    if save or stub:
        _save_transferred_record(record)
    if stub:
        return registry.get(uid=record.uid)
    return None
