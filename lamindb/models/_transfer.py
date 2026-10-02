from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any

import lamindb_setup as ln_setup
from django.core.exceptions import FieldDoesNotExist
from django.db import ProgrammingError
from django.db.models import Model
from django.db.models import QuerySet as DjangoQuerySet
from lamin_utils import logger
from lamindb_setup._connect_instance import get_owner_name_from_identifier
from lamindb_setup.errors import NoReadAccess

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
            field = REGISTRY_UNIQUE_FIELD.get(fk, "uid")
            if field == "uid":
                pre_existing_fk_record_default = _cached_or_load(
                    fk_record, transfer_logs
                )
            else:
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
    # The transfer run is the lineage. An ambient ln.track() run, when present,
    # is only the parent (initiated_by_run). `lamin io sync` has no such parent.
    initiated_by_run = context.run
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


@dataclass
class AnnotationGap:
    """A feature link whose feature or value is not visible to this account.

    ``feature_name`` is set only when the feature row itself can be read.
    ``hidden_value_ids`` are ids already stored on the link table.
    """

    feature_id: int | None
    feature_name: str | None
    feature_uid: str | None
    value_model: str | None
    hidden_value_ids: list[int]
    n_total: int
    n_hidden: int


def _annotation_value_links(record):
    """Yield ``(accessor, value_field)`` for feature-value link tables.

    Link tables have no space of their own, so their rows stay visible when the
    value (or the feature) is in a space this account cannot read.
    """
    for rel in record._meta.related_objects:
        accessor = rel.get_accessor_name()
        if not accessor or not str(accessor).startswith("values_"):
            continue
        model = rel.related_model
        try:
            value_field = model._meta.get_field("value")
        except FieldDoesNotExist:
            # Run.values_artifact stores the target as `artifact`, not `value`.
            continue
        yield accessor, value_field


def _blocked_annotation_message(record, gaps: list[AnnotationGap]) -> str:
    lines = [
        f"Cannot transfer annotations of {record.__class__.__name__} {record.uid!r}."
    ]
    for gap in gaps:
        hidden_ids = ", ".join(str(value_id) for value_id in gap.hidden_value_ids)
        id_label = "id" if len(gap.hidden_value_ids) == 1 else "ids"
        partial = ""
        if gap.n_hidden != gap.n_total:
            partial = f" ({gap.n_hidden} of {gap.n_total} values)"
        lines.append(
            f"Feature {gap.feature_name!r} (uid={gap.feature_uid}) links a "
            f"{gap.value_model} ({id_label}={hidden_ids}){partial} that this "
            "account cannot read."
        )
    lines.append("The annotation set would be incomplete.")
    lines.append('Pass transfer="sqlrecord" to sync the object without annotations.')
    return "\n".join(lines)


_NOT_CACHED = object()


def _pop_cached_linked_values(transfer_logs: dict, record):
    cache = transfer_logs.get("_linked_values")
    if not isinstance(cache, dict):
        return None
    return cache.pop(getattr(record, "uid", None), None)


def _cached_or_load(record, transfer_logs: dict):
    """Target row for ``record.uid``, looked up once per registry batch.

    ``None`` means the uid was looked up and is absent. A missing cache entry
    loads that one uid and stores the result.
    """
    uid = getattr(record, "uid", None)
    if uid is None:
        return None
    resolved = transfer_logs.setdefault("_resolved", {})
    key = (record.__class__.__name__, uid)
    if key in resolved:
        return resolved[key]
    found = record.__class__.objects.filter(uid=uid).one_or_none()
    resolved[key] = found
    return found


def _remember_target(record, transfer_logs: dict) -> None:
    uid = getattr(record, "uid", None)
    if uid is None:
        return
    transfer_logs.setdefault("_resolved", {})[(record.__class__.__name__, uid)] = record


def resolve_records(records, transfer_logs: dict) -> None:
    """One ``uid__in`` lookup per registry for records that are not cached yet."""
    bucket: dict = {}
    for record in records:
        uid = getattr(record, "uid", None)
        if uid is None:
            continue
        bucket.setdefault(record.__class__, {})[uid] = record
    _resolve_present(bucket, transfer_logs)


def _resolve_present(bucket: dict, transfer_logs: dict) -> None:
    resolved = transfer_logs.setdefault("_resolved", {})
    for model, by_uid in bucket.items():
        unknown = [uid for uid in by_uid if (model.__name__, uid) not in resolved]
        if not unknown:
            continue
        found = {row.uid: row for row in model.objects.filter(uid__in=unknown)}
        for uid in unknown:
            resolved[(model.__name__, uid)] = found.get(uid)


def _put_entity(bucket: dict, record) -> None:
    uid = getattr(record, "uid", None)
    if uid is None:
        return
    if not isinstance(record, (SQLRecord, BaseSQLRecord)):
        return
    bucket.setdefault(record.__class__, {})[uid] = record


def _collect_entities(linked_values, bucket: dict) -> None:
    for feature, value in linked_values:
        _put_entity(bucket, feature)
        values = value if isinstance(value, list) else [value]
        for item in values:
            _put_entity(bucket, item)


def _read_link_rows(record) -> tuple[list[dict], dict[str, type[Model]]]:
    """All feature-value links for one record, in one query.

    Relational rows carry ``value_pk``. JSON rows carry ``json_value``.
    A ``value_pk`` that does not resolve is an unreadable annotation.
    """
    from django.db.models import BigIntegerField, CharField, F, JSONField, Value

    model_by_name: dict[str, type[Model]] = {}
    parts = []
    for accessor, value_field in _annotation_value_links(record):
        relational = bool(getattr(value_field, "is_relation", False))
        model_name = value_field.related_model.__name__ if relational else ""
        if relational:
            model_by_name[model_name] = value_field.related_model
        qs = getattr(record, accessor).order_by()
        if relational:
            qs = qs.annotate(
                link_id=F("id"),
                value_pk=F("value_id"),
                json_value=Value(None, output_field=JSONField()),
                value_model=Value(model_name, output_field=CharField(max_length=64)),
            )
        else:
            qs = qs.annotate(
                link_id=F("id"),
                value_pk=Value(None, output_field=BigIntegerField()),
                json_value=F("value"),
                value_model=Value("", output_field=CharField(max_length=64)),
            )
        parts.append(
            qs.values("link_id", "feature_id", "value_pk", "json_value", "value_model")
        )
    if not parts:
        return [], model_by_name
    rows = parts[0] if len(parts) == 1 else parts[0].union(*parts[1:], all=True)
    return list(rows), model_by_name


def _assemble_linked_values(record, rows: list[dict], features: dict, values: dict):
    """Group link rows into ``(feature, value)`` pairs.

    Several categorical links for one feature are one list. A JSON list stays
    one value because it is stored as a single JSON cell.
    """
    rows = sorted(rows, key=lambda row: row["link_id"] or 0)
    grouped: dict[int, list] = {}
    feature_for: dict[int, Any] = {}
    totals: dict[int, int] = {}
    gap_slots: dict[int, dict] = {}
    for row in rows:
        feature_id = row["feature_id"]
        feature = features.get(feature_id)
        if feature is None:
            raise NoReadAccess(
                _blocked_annotation_message(
                    record,
                    [
                        AnnotationGap(
                            feature_id=feature_id,
                            feature_name=None,
                            feature_uid=None,
                            value_model=row["value_model"] or None,
                            hidden_value_ids=[row["value_pk"]]
                            if row["value_pk"] is not None
                            else [],
                            n_total=1,
                            n_hidden=1,
                        )
                    ],
                )
            )
        feature_for[feature_id] = feature
        model_name = row["value_model"] or ""
        if not model_name:
            grouped.setdefault(feature_id, []).append(row["json_value"])
            continue
        totals[feature_id] = totals.get(feature_id, 0) + 1
        value = values.get((model_name, row["value_pk"]))
        if value is None:
            slot = gap_slots.setdefault(
                feature_id,
                {"feature": feature, "model": model_name, "hidden_ids": []},
            )
            slot["hidden_ids"].append(row["value_pk"])
            continue
        grouped.setdefault(feature_id, []).append(value)
    if gap_slots:
        gaps = []
        for feature_id, slot in gap_slots.items():
            feature = slot["feature"]
            hidden_ids = list(slot["hidden_ids"])
            gaps.append(
                AnnotationGap(
                    feature_id=feature.id,
                    feature_name=feature.name,
                    feature_uid=feature.uid,
                    value_model=slot["model"],
                    hidden_value_ids=hidden_ids,
                    n_total=totals[feature_id],
                    n_hidden=len(hidden_ids),
                )
            )
        raise NoReadAccess(_blocked_annotation_message(record, gaps))
    return [
        (feature_for[feature_id], vals[0] if len(vals) == 1 else vals)
        for feature_id, vals in grouped.items()
    ]


def _hydrate_link_rows(record, rows: list[dict], model_by_name: dict[str, type[Model]]):
    from .feature import Feature

    source_db = record._state.db
    feature_ids = {row["feature_id"] for row in rows if row["feature_id"] is not None}
    features = {}
    if feature_ids:
        features = {
            feature.id: feature
            for feature in Feature.objects.using(source_db).filter(id__in=feature_ids)
        }
    value_ids: dict[str, set] = {}
    for row in rows:
        if row["value_model"] and row["value_pk"] is not None:
            value_ids.setdefault(row["value_model"], set()).add(row["value_pk"])
    values = {}
    for model_name, ids in value_ids.items():
        model = model_by_name[model_name]
        for obj in model.objects.using(source_db).filter(id__in=ids):
            values[(model_name, obj.id)] = obj
    return _assemble_linked_values(record, rows, features, values)


def _linked_feature_values(record) -> list[tuple[Any, Any]]:
    """Feature values from link rows, keyed by the feature row rather than its name.

    Names are not unique. Several categorical links for one feature are one list.
    A JSON list stays one value because it is stored as a single JSON cell.

    Raises ``NoReadAccess`` when a link points at a value this account cannot
    read. Skipping it would transfer a partial annotation set.
    """
    rows, model_by_name = _read_link_rows(record)
    if not rows:
        return []
    return _hydrate_link_rows(record, rows, model_by_name)


def prime_annotation_transfer(hosts: list, transfer_logs: dict) -> None:
    """Read every host's links and resolve those uids on the target once.

    One link query per host. Feature rows and value rows are loaded with one
    ``id__in`` per registry on the source, then one ``uid__in`` per registry
    on the target.
    """
    cache = transfer_logs.setdefault("_linked_values", {})
    pending = []
    for host in hosts:
        if host.__class__.__name__ not in {"Record", "Run"}:
            continue
        if host.uid in cache:
            continue
        rows, model_by_name = _read_link_rows(host)
        pending.append((host, rows, model_by_name))
    if not pending:
        return
    by_db: dict = {}
    for host, rows, model_by_name in pending:
        by_db.setdefault(host._state.db, []).append((host, rows, model_by_name))
    bucket: dict = {}
    for source_db, items in by_db.items():
        _hydrate_hosts(source_db, items, cache, bucket)
    _resolve_present(bucket, transfer_logs)


def _hydrate_hosts(source_db, items, cache: dict, bucket: dict) -> None:
    from .feature import Feature

    feature_ids: set = set()
    value_ids: dict[str, set] = {}
    model_by_name: dict[str, type[Model]] = {}
    for _host, rows, names in items:
        model_by_name.update(names)
        for row in rows:
            if row["feature_id"] is not None:
                feature_ids.add(row["feature_id"])
            if row["value_model"] and row["value_pk"] is not None:
                value_ids.setdefault(row["value_model"], set()).add(row["value_pk"])
    features = {}
    if feature_ids:
        features = {
            feature.id: feature
            for feature in Feature.objects.using(source_db).filter(id__in=feature_ids)
        }
    values = {}
    for model_name, ids in value_ids.items():
        model = model_by_name[model_name]
        for obj in model.objects.using(source_db).filter(id__in=ids):
            values[(model_name, obj.id)] = obj
    for host, rows, _names in items:
        linked = _assemble_linked_values(host, rows, features, values)
        cache[host.uid] = linked
        _collect_entities(linked, bucket)


def _clear_annotation_links(record, transfer_logs: dict) -> None:
    """Delete this record's annotation links, one query per link table.

    A row inserted earlier in this sync has no links to replace.
    """
    inserted = transfer_logs.get("_inserted")
    if isinstance(inserted, set) and record.uid in inserted:
        return
    if record.pk is None:
        return
    for rel in record._meta.related_objects:
        accessor = rel.get_accessor_name()
        if not accessor or not str(accessor).startswith("values_"):
            continue
        try:
            rel.related_model._meta.get_field("value")
        except FieldDoesNotExist:
            continue
        rel.related_model.objects.filter(**{rel.field.name: record.pk}).delete()


def _depth_descendants(record, depth: int) -> list:
    from .sqlrecord import _typed_children

    if depth <= 0:
        return []
    children = _typed_children(record)
    found = list(children)
    for child in children:
        found.extend(_depth_descendants(child, depth - 1))
    return found


def log_transferred_record(
    record, transfer_logs: dict, mapped_before: int, transferred_before: int
) -> None:
    name = (
        getattr(record, "name", None)
        or getattr(record, "key", None)
        or getattr(record, "uid", None)
    )
    n_new = len(transfer_logs["transferred"]) - transferred_before
    n_have = len(transfer_logs["mapped"]) - mapped_before
    logger.important(
        f"{type(record).__name__} {name}: {n_new} transferred, {n_have} already on target"
    )


def transfer_record_feature_values(
    record_on_default, source_db, source_pk, using, transfer_logs
):
    from copy import copy

    from .feature import Feature, parse_dtype

    linked_values = _pop_cached_linked_values(transfer_logs, record_on_default)
    if linked_values is None:
        source = record_on_default.__class__.objects.using(source_db).get(pk=source_pk)
        linked_values = _linked_feature_values(source)
    if not linked_values:
        return

    def _transfer_entity(value, feature=None):
        known = transfer_logs.get("_resolved", {}).get(
            (type(value).__name__, getattr(value, "uid", None)), _NOT_CACHED
        )
        if known is not _NOT_CACHED and known is not None:
            transfer_logs["mapped"].append(f"{type(value).__name__}(uid='{value.uid}')")
            if type(value).__name__ == "User":
                return getattr(known, _user_annotation_field(feature))
            return known
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
        return value.save(
            transfer="annotations",
            _transfer_logs=transfer_logs,
            _transfer_summarize=False,
        )

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
    _clear_annotation_links(record_on_default, transfer_logs)
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
    # Read every annotation target before writing the row. Link tables are
    # visible even when the value sits in a space this account cannot read,
    # and saving first would leave a record whose annotations are not the source.
    if (
        transfer_annotations
        and not stub
        and record.__class__.__name__ in {"Record", "Run"}
    ):
        cache = transfer_logs.setdefault("_linked_values", {})
        if record.uid not in cache:
            cache[record.uid] = _linked_feature_values(record)
    # Dtype text is not a foreign key. Follow it even when this feature row
    # is already on the target, so a re-transfer picks up schema__uid refs.
    if record.__class__.__name__ == "Feature":
        from .feature import transfer_feature_dtypes

        transfer_feature_dtypes(record, using, transfer_logs=transfer_logs)
    registry = record.__class__
    logger.debug(f"transferring {registry.__name__} record {record.uid} to default db")
    record_on_default = _cached_or_load(record, transfer_logs)
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

        record = copy(record)
    else:
        transfer_logs["transferred"].append(record_str)
        # The caller may insert this row itself. It has no annotation links yet.
        transfer_logs.setdefault("_inserted", set()).add(record.uid)

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
        _remember_target(record_on_default, transfer_logs)
        return record_on_default
    record.id = None
    if save or stub:
        _save_transferred_record(record)
        transfer_logs.setdefault("_inserted", set()).add(record.uid)
        saved_row = registry.get(uid=record.uid) if stub else record
        _remember_target(saved_row, transfer_logs)
        if stub:
            return saved_row
    return None


def _registry_class_name(registry: str) -> str:
    if registry == "ulabel":
        return "ULabel"
    if not registry or not registry.replace("_", "").isalnum():
        raise ValueError(f"Unknown registry {registry!r}.")
    return "".join(part.capitalize() for part in registry.split("_"))


def sync_objects_from_database(
    registry: str,
    uids: str | list[str],
    *,
    source: str,
    depth: int = 0,
    transfer: str | None = None,
) -> list[SQLRecord]:
    """Sync SQLRecord objects from a source database into the default database.

    This is a high-level function used in the CLI: `lamin io sync`.

    One sync walks a single graph. Shared nodes are copied once. Link rows are
    replaced on the target, not copied by uid.

    .. code-block:: mermaid

       flowchart TD
         roots["Requested records and depth children"] --> row["Copy the row once; fill Record and ULabel"]
         row --> shared["Schema, features, dtype types: once per uid"]
         row --> lookup["Annotation values: one uid lookup per registry"]
         lookup --> present["Already on target: use that row"]
         lookup --> stub["Missing Record or ULabel: stub"]
         lookup --> once["Missing other registry: save once"]
         row --> links["Replace link rows; do not copy them by uid"]

    Most of the time, you will just `.save()` on an object from another database::

        import lamindb as ln
        db = ln.DB("laminlabs/lamindata")
        record = db.Record.get(uid="gL3TbX2qZQmCwTAU")
        record.save(transfer="sqlrecord")

    Guide: {doc}`transfer`

    Args:
        registry: Registry name, for example `artifact` or `record`.
        uids: One uid or several uids on the source database.
        source: Source instance slug, for example `laminlabs/lamindata`.
        depth: How many levels of records under a type to transfer.
            `0` transfers only the given objects, plus related objects selected
            by `transfer`. A positive integer also transfers that many levels
            of records whose type chain starts at each object. Only `record`,
            `feature`, `schema`, `project`, `ulabel`, and `reference` accept
            `depth > 0`.
        transfer: `sqlrecord`, `notes`, or `annotations`.
            Omit it to use the registry default.
    """
    from .db import DB

    if type(depth) is not int or depth < 0:
        raise ValueError("depth must be an int >= 0.")
    if isinstance(uids, str):
        uid_list = [uids]
    elif isinstance(uids, list):
        uid_list = list(uids)
    else:
        raise TypeError("uids must be a str or list[str].")
    if not uid_list:
        raise ValueError("uids is required and must contain at least one uid.")
    model_name = _registry_class_name(registry)
    queryset = getattr(DB(source), model_name)
    saved: list[SQLRecord] = []
    for uid in uid_list:
        record = queryset.get(uid)
        kwargs: dict[str, Any] = {"depth": depth}
        if transfer is not None:
            kwargs["transfer"] = transfer
        saved.append(record.save(**kwargs))
    return saved
