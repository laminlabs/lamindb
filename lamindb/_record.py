from __future__ import annotations

import builtins
from typing import TYPE_CHECKING, List, NamedTuple

import dj_database_url
import lamindb_setup as ln_setup
from django.db import connections
from django.db.models import IntegerField, Manager, Q, QuerySet, Value
from lamin_utils import logger
from lamin_utils._lookup import Lookup
from lamindb_setup._connect_instance import get_owner_name_from_identifier
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core._hub_core import connect_instance
from lnschema_core.models import IsVersioned, Record

from lamindb._utils import attach_func_to_class_method
from lamindb.core._settings import settings

from ._from_values import get_or_create_records

if TYPE_CHECKING:
    import pandas as pd
    from lnschema_core.types import ListLike, StrField


IPYTHON = getattr(builtins, "__IPYTHON__", False)


def init_self_from_db(self: Record, existing_record: Record):
    new_args = [
        getattr(existing_record, field.attname) for field in self._meta.concrete_fields
    ]
    super(self.__class__, self).__init__(*new_args)
    self._state.adding = False  # mimic from_db
    self._state.db = "default"


def validate_required_fields(orm: Record, kwargs):
    required_fields = {
        k.name for k in orm._meta.fields if not k.null and k.default is None
    }
    required_fields_not_passed = {k: None for k in required_fields if k not in kwargs}
    kwargs.update(required_fields_not_passed)
    missing_fields = [
        k for k, v in kwargs.items() if v is None and k in required_fields
    ]
    if missing_fields:
        raise TypeError(f"{missing_fields} are required.")


def suggest_records_with_similar_names(record: Record, kwargs) -> bool:
    """Returns True if found exact match, otherwise False.

    Logs similar matches if found.
    """
    if kwargs.get("name") is None:
        return False
    queryset = _search(
        record.__class__, kwargs["name"], field="name", truncate_words=True, limit=20
    )
    if not queryset.exists():  # empty queryset
        return False
    for alternative_record in queryset:
        if alternative_record.name == kwargs["name"]:
            return True
    s, it, nots = ("", "it", "s") if len(queryset) == 1 else ("s", "one of them", "")
    msg = f"record{s} with similar name{s} exist{nots}! did you mean to load {it}?"
    if IPYTHON:
        from IPython.display import display

        logger.warning(f"{msg}")
        if settings._verbosity_int >= 1:
            display(queryset.df())
    else:
        logger.warning(f"{msg}\n{queryset}")
    return False


def __init__(orm: Record, *args, **kwargs):
    if not args:
        validate_required_fields(orm, kwargs)

        # do not search for names if an id is passed; this is important
        # e.g. when synching ids from the notebook store to lamindb
        has_consciously_provided_uid = False
        if "_has_consciously_provided_uid" in kwargs:
            has_consciously_provided_uid = kwargs.pop("_has_consciously_provided_uid")
        if settings.creation.search_names and not has_consciously_provided_uid:
            match = suggest_records_with_similar_names(orm, kwargs)
            if match:
                if "version" in kwargs:
                    version_comment = " and version"
                    existing_record = orm.__class__.filter(
                        name=kwargs["name"], version=kwargs["version"]
                    ).one_or_none()
                else:
                    version_comment = ""
                    existing_record = orm.__class__.filter(name=kwargs["name"]).one()
                if existing_record is not None:
                    logger.important(
                        f"returning existing {orm.__class__.__name__} record with same"
                        f" name{version_comment}: '{kwargs['name']}'"
                    )
                    init_self_from_db(orm, existing_record)
                    return None
        super(Record, orm).__init__(**kwargs)
    elif len(args) != len(orm._meta.concrete_fields):
        raise ValueError("please provide keyword arguments, not plain arguments")
    else:
        # object is loaded from DB (**kwargs could be omitted below, I believe)
        super(Record, orm).__init__(*args, **kwargs)


@classmethod  # type:ignore
@doc_args(Record.filter.__doc__)
def filter(cls, **expressions) -> QuerySet:
    """{}"""  # noqa: D415
    from lamindb._filter import filter

    return filter(cls, **expressions)


@classmethod  # type:ignore
@doc_args(Record.get.__doc__)
def get(cls, idlike: int | str) -> Record:
    """{}"""  # noqa: D415
    from lamindb._filter import filter

    if isinstance(idlike, int):
        return filter(cls, id=idlike).one()
    else:
        qs = filter(cls, uid__startswith=idlike)
        if issubclass(cls, IsVersioned):
            return qs.latest_version().one()
        else:
            return qs.one()


@classmethod  # type:ignore
@doc_args(Record.df.__doc__)
def df(
    cls, include: str | list[str] | None = None, join: str = "inner"
) -> pd.DataFrame:
    """{}"""  # noqa: D415
    from lamindb._filter import filter

    query_set = filter(cls)
    if hasattr(cls, "updated_at"):
        query_set = query_set.order_by("-updated_at")
    return query_set.df(include=include, join=join)


# from_values doesn't apply for QuerySet or Manager
@classmethod  # type:ignore
@doc_args(Record.from_values.__doc__)
def from_values(
    cls,
    values: ListLike,
    field: StrField | None = None,
    create: bool = False,
    organism: Record | str | None = None,
    public_source: Record | None = None,
    mute: bool = False,
) -> list[Record]:
    """{}"""  # noqa: D415
    from_public = True if cls.__module__.startswith("lnschema_bionty.") else False
    field_str = get_default_str_field(cls, field=field)
    return get_or_create_records(
        iterable=values,
        field=getattr(cls, field_str),
        create=create,
        from_public=from_public,
        organism=organism,
        public_source=public_source,
        mute=mute,
    )


def _search(
    cls,
    string: str,
    *,
    field: StrField | list[StrField] | None = None,
    limit: int | None = 20,
    case_sensitive: bool = False,
    using_key: str | None = None,
    truncate_words: bool = False,
) -> QuerySet:
    input_queryset = _queryset(cls, using_key=using_key)
    orm = input_queryset.model
    if field is None:
        fields = [
            field.name
            for field in orm._meta.fields
            if field.get_internal_type() in {"CharField", "TextField"}
        ]
    else:
        if not isinstance(field, list):
            fields_input = [field]
        else:
            fields_input = field
        fields = []
        for field in fields_input:
            if not isinstance(field, str):
                try:
                    fields.append(field.field.name)
                except AttributeError as error:
                    raise TypeError(
                        "Please pass a Record string field, e.g., `CellType.name`!"
                    ) from error
            else:
                fields.append(field)

    # decompose search string
    def truncate_word(word) -> str:
        if len(word) > 5:
            n_80_pct = int(len(word) * 0.8)
            return word[:n_80_pct]
        elif len(word) > 3:
            return word[:3]
        else:
            return word

    decomposed_string = str(string).split()
    # add the entire string back
    decomposed_string += [string]
    for word in decomposed_string:
        # will not search against words with 3 or fewer characters
        if len(word) <= 3:
            decomposed_string.remove(word)
    if truncate_words:
        decomposed_string = [truncate_word(word) for word in decomposed_string]
    # construct the query
    expression = Q()
    case_sensitive_i = "" if case_sensitive else "i"
    for field in fields:
        for word in decomposed_string:
            query = {f"{field}__{case_sensitive_i}contains": word}
            expression |= Q(**query)
    output_queryset = input_queryset.filter(expression)
    # ensure exact matches are at the top
    narrow_expression = Q()
    for field in fields:
        query = {f"{field}__{case_sensitive_i}contains": string}
        narrow_expression |= Q(**query)
    refined_output_queryset = output_queryset.filter(narrow_expression).curate(
        ordering=Value(1, output_field=IntegerField())
    )
    remaining_output_queryset = output_queryset.exclude(narrow_expression).curate(
        ordering=Value(2, output_field=IntegerField())
    )
    combined_queryset = refined_output_queryset.union(
        remaining_output_queryset
    ).order_by("ordering")[:limit]
    return combined_queryset


@classmethod  # type: ignore
@doc_args(Record.search.__doc__)
def search(
    cls,
    string: str,
    *,
    field: StrField | None = None,
    limit: int | None = 20,
    case_sensitive: bool = False,
) -> QuerySet:
    """{}"""  # noqa: D415
    return _search(
        cls=cls,
        string=string,
        field=field,
        limit=limit,
        case_sensitive=case_sensitive,
    )


def _lookup(
    cls,
    field: StrField | None = None,
    return_field: StrField | None = None,
    using_key: str | None = None,
) -> NamedTuple:
    """{}"""  # noqa: D415
    queryset = _queryset(cls, using_key=using_key)
    field = get_default_str_field(orm=queryset.model, field=field)

    return Lookup(
        records=queryset,
        values=[i.get(field) for i in queryset.values()],
        tuple_name=cls.__class__.__name__,
        prefix="ln",
    ).lookup(
        return_field=(
            get_default_str_field(orm=queryset.model, field=return_field)
            if return_field is not None
            else None
        )
    )


@classmethod  # type: ignore
@doc_args(Record.lookup.__doc__)
def lookup(
    cls,
    field: StrField | None = None,
    return_field: StrField | None = None,
) -> NamedTuple:
    """{}"""  # noqa: D415
    return _lookup(cls=cls, field=field, return_field=return_field)


def get_default_str_field(
    orm: Record | QuerySet | Manager,
    *,
    field: str | StrField | None = None,
) -> str:
    """Get the 1st char or text field from the orm."""
    if isinstance(orm, (QuerySet, Manager)):
        orm = orm.model
    model_field_names = [i.name for i in orm._meta.fields]

    # set default field
    if field is None:
        if orm._meta.model.__name__ == "Run":
            field = orm._meta.get_field("created_at")
        elif orm._meta.model.__name__ == "User":
            field = orm._meta.get_field("handle")
        elif "name" in model_field_names:
            # by default use the name field
            field = orm._meta.get_field("name")
        else:
            # first char or text field that doesn't contain "id"
            for i in orm._meta.fields:
                if "id" in i.name:
                    continue
                if i.get_internal_type() in {"CharField", "TextField"}:
                    field = i
                    break

        # no default field can be found
        if field is None:
            raise ValueError(
                "please pass a Record string field, e.g., `CellType.name`!"
            )
        else:
            field = field.name  # type:ignore
    if not isinstance(field, str):
        try:
            field = field.field.name
        except AttributeError:
            raise TypeError(
                "please pass a Record string field, e.g., `CellType.name`!"
            ) from None

    return field


def _queryset(cls: Record | QuerySet | Manager, using_key: str) -> QuerySet:
    if isinstance(cls, (QuerySet, Manager)):
        return cls.all()
    elif using_key is None:
        return cls.objects.all()
    else:
        # using must be called on cls, otherwise the connection isn't found
        return cls.using(using_key).all()


def add_db_connection(db: str, using: str):
    db_config = dj_database_url.config(
        default=db, conn_max_age=600, conn_health_checks=True
    )
    db_config["TIME_ZONE"] = "UTC"
    db_config["OPTIONS"] = {}
    db_config["AUTOCOMMIT"] = True
    connections.settings[using] = db_config


@classmethod  # type: ignore
@doc_args(Record.using.__doc__)
def using(
    cls,
    instance: str,
) -> QuerySet:
    """{}"""  # noqa: D415
    from lamindb_setup._connect_instance import (
        load_instance_settings,
        update_db_using_local,
    )
    from lamindb_setup.core._settings_store import instance_settings_file

    owner, name = get_owner_name_from_identifier(instance)
    settings_file = instance_settings_file(name, owner)
    if not settings_file.exists():
        load_result = connect_instance(owner=owner, name=name)
        if isinstance(load_result, str):
            raise RuntimeError(
                f"Failed to load instance {instance}, please check your permission!"
            )
        instance_result, _ = load_result
        settings_file = instance_settings_file(name, owner)
        db = update_db_using_local(instance_result, settings_file)
    else:
        isettings = load_instance_settings(settings_file)
        db = isettings.db
    add_db_connection(db, instance)
    return QuerySet(model=cls, using=instance)


REGISTRY_UNIQUE_FIELD = {
    "storage": "root",
    "feature": "name",
    "ulabel": "name",
}


def update_fk_to_default_db(
    records: Record | list[Record] | QuerySet,
    fk: str,
    using_key: str | None,
):
    record = records[0] if isinstance(records, (List, QuerySet)) else records
    if hasattr(record, f"{fk}_id") and getattr(record, f"{fk}_id") is not None:
        fk_record = getattr(record, fk)
        field = REGISTRY_UNIQUE_FIELD.get(fk, "uid")
        fk_record_default = fk_record.__class__.filter(
            **{field: getattr(fk_record, field)}
        ).one_or_none()
        if fk_record_default is None:
            from copy import copy

            fk_record_default = copy(fk_record)
            transfer_to_default_db(fk_record_default, using_key, save=True)
        if isinstance(records, (List, QuerySet)):
            for r in records:
                setattr(r, f"{fk}", None)
                setattr(r, f"{fk}_id", fk_record_default.id)
        else:
            setattr(records, f"{fk}", None)
            setattr(records, f"{fk}_id", fk_record_default.id)


FKBULK = [
    "organism",
    "public_source",
    "latest_report",  # Transform
    "source_code",  # Transform
    "report",  # Run
]


def transfer_fk_to_default_db_bulk(records: list | QuerySet, using_key: str | None):
    for fk in FKBULK:
        update_fk_to_default_db(records, fk, using_key)


def transfer_to_default_db(
    record: Record,
    using_key: str | None,
    save: bool = False,
    mute: bool = False,
    transfer_fk: bool = True,
) -> Record | None:
    db = record._state.db
    if db is not None and db != "default" and using_key is None:
        registry = record.__class__
        record_on_default = registry.objects.filter(uid=record.uid).one_or_none()
        if record_on_default is not None:
            logger.important(
                f"returning existing {record.__class__.__name__}(uid='{record.uid}') on default database"
            )
            return record_on_default
        if not mute:
            logger.hint(f"saving from instance {db} to default instance: {record}")
        from lamindb.core._data import WARNING_RUN_TRANSFORM
        from lamindb.core._run_context import run_context

        if hasattr(record, "created_by_id"):
            # this line is needed to point created_by to default db
            record.created_by = None
            record.created_by_id = ln_setup.settings.user.id
        if hasattr(record, "run_id"):
            record.run = None
            if run_context.run is not None:
                record.run_id = run_context.run.id
            else:
                if not settings.creation.artifact_silence_missing_run_warning:
                    logger.warning(WARNING_RUN_TRANSFORM)
                record.run_id = None
        if hasattr(record, "transform_id") and record._meta.model_name != "run":
            record.transform = None
            if run_context.transform is not None:
                record.transform_id = run_context.transform.id
            else:
                record.transform_id = None
        # transfer other foreign key fields
        fk_fields = [
            i.name
            for i in record._meta.fields
            if i.get_internal_type() == "ForeignKey"
            if i.name not in {"created_by", "run", "transform"}
        ]
        if not transfer_fk:
            # don't transfer fk fields that are already bulk transferred
            fk_fields = [fk for fk in fk_fields if fk not in FKBULK]
        for fk in fk_fields:
            update_fk_to_default_db(record, fk, using_key)
        record.id = None
        record._state.db = "default"
        if save:
            record.save()
    return None


# docstring handled through attach_func_to_class_method
def save(self, *args, **kwargs) -> Record:
    using_key = None
    if "using" in kwargs:
        using_key = kwargs["using"]
    db = self._state.db
    pk_on_db = self.pk
    artifacts: list = []
    if self.__class__.__name__ == "Collection" and self.id is not None:
        # when creating a new collection without being able to access artifacts
        artifacts = self.artifacts.list()
    # transfer of the record to the default db with fk fields
    result = transfer_to_default_db(self, using_key)
    if result is not None:
        init_self_from_db(self, result)
    else:
        # here, we can't use the parents argument
        # parents are not saved for the self record
        save_kwargs = kwargs.copy()
        if "parents" in save_kwargs:
            save_kwargs.pop("parents")
        super(Record, self).save(*args, **save_kwargs)
    # perform transfer of many-to-many fields
    # only supported for Artifact and Collection records
    if db is not None and db != "default" and using_key is None:
        if self.__class__.__name__ == "Collection":
            if len(artifacts) > 0:
                logger.info("transfer artifacts")
                for artifact in artifacts:
                    artifact.save()
                self.unordered_artifacts.add(*artifacts)
        if hasattr(self, "labels"):
            from copy import copy

            from lnschema_core.models import FeatureManager

            # here we go back to original record on the source database
            self_on_db = copy(self)
            self_on_db._state.db = db
            self_on_db.pk = pk_on_db  # manually set the primary key
            self_on_db.features = FeatureManager(self_on_db)
            # by default, transfer parents of the labels to maintain ontological hierarchy
            try:
                import bionty as bt

                parents = kwargs.get("parents", bt.settings.auto_save_parents)
            except ImportError:
                parents = kwargs.get("parents", True)
            add_from_kwargs = {"parents": parents}
            self.features._add_from(self_on_db, **add_from_kwargs)
            self.labels.add_from(self_on_db, **add_from_kwargs)
    return self


METHOD_NAMES = [
    "__init__",
    "filter",
    "get",
    "df",
    "search",
    "lookup",
    "save",
    "from_values",
    "using",
]

if ln_setup._TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(Record, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Record, globals())
