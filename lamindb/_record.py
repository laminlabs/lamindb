from __future__ import annotations

import builtins
from typing import TYPE_CHECKING, NamedTuple

import dj_database_url
import lamindb_setup as ln_setup
from django.core.exceptions import FieldDoesNotExist
from django.db import connections, transaction
from django.db.models import IntegerField, Manager, Q, QuerySet, Value
from lamin_utils import colors, logger
from lamin_utils._lookup import Lookup
from lamindb_setup._connect_instance import (
    get_owner_name_from_identifier,
    load_instance_settings,
    update_db_using_local,
)
from lamindb_setup.core._docs import doc_args
from lamindb_setup.core._hub_core import connect_instance_hub
from lamindb_setup.core._settings_store import instance_settings_file
from lnschema_core.models import Artifact, Feature, IsVersioned, Record, Run, Transform

from ._utils import attach_func_to_class_method
from .core._settings import settings
from .core.exceptions import RecordNameChangeIntegrityError, ValidationError

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


def update_attributes(record: Record, attributes: dict[str, str]):
    for key, value in attributes.items():
        if getattr(record, key) != value:
            logger.warning(f"updated {key} from {getattr(record, key)} to {value}")
            setattr(record, key, value)


def validate_required_fields(record: Record, kwargs):
    required_fields = {
        k.name for k in record._meta.fields if not k.null and k.default is None
    }
    required_fields_not_passed = {k: None for k in required_fields if k not in kwargs}
    kwargs.update(required_fields_not_passed)
    missing_fields = [
        k for k, v in kwargs.items() if v is None and k in required_fields
    ]
    if missing_fields:
        raise TypeError(f"{missing_fields} are required.")
    try:
        uid_max_length = record.__class__._meta.get_field(
            "uid"
        ).max_length  # triggers FieldDoesNotExist
        if len(kwargs["uid"]) != uid_max_length:  # triggers KeyError
            raise ValidationError(
                f'`uid` must be exactly {uid_max_length} characters long, got {len(kwargs["uid"])}.'
            )
    except (FieldDoesNotExist, KeyError):
        pass


def suggest_records_with_similar_names(record: Record, name_field: str, kwargs) -> bool:
    """Returns True if found exact match, otherwise False.

    Logs similar matches if found.
    """
    if kwargs.get(name_field) is None or not isinstance(kwargs.get(name_field), str):
        return False
    queryset = _search(
        record.__class__,
        kwargs[name_field],
        field=name_field,
        truncate_words=True,
        limit=3,
    )
    if not queryset.exists():  # empty queryset
        return False
    for alternative_record in queryset:
        if getattr(alternative_record, name_field) == kwargs[name_field]:
            return True
    s, it, nots = ("", "it", "s") if len(queryset) == 1 else ("s", "one of them", "")
    msg = f"record{s} with similar {name_field}{s} exist{nots}! did you mean to load {it}?"
    if IPYTHON:
        from IPython.display import display

        logger.warning(f"{msg}")
        if settings._verbosity_int >= 1:
            display(queryset.df())
    else:
        logger.warning(f"{msg}\n{queryset}")
    return False


def __init__(record: Record, *args, **kwargs):
    if not args:
        validate_required_fields(record, kwargs)

        # do not search for names if an id is passed; this is important
        # e.g. when synching ids from the notebook store to lamindb
        has_consciously_provided_uid = False
        if "_has_consciously_provided_uid" in kwargs:
            has_consciously_provided_uid = kwargs.pop("_has_consciously_provided_uid")
        if settings.creation.search_names and not has_consciously_provided_uid:
            name_field = (
                "name" if not hasattr(record, "_name_field") else record._name_field
            )
            match = suggest_records_with_similar_names(record, name_field, kwargs)
            if match:
                if "version" in kwargs:
                    if kwargs["version"] is not None:
                        version_comment = " and version"
                        existing_record = record.__class__.filter(
                            **{
                                name_field: kwargs[name_field],
                                "version": kwargs["version"],
                            }
                        ).one_or_none()
                    else:
                        # for a versioned record, an exact name match is not a
                        # criterion for retrieving a record in case `version`
                        # isn't passed - we'd always pull out many records with exactly the
                        # same name
                        existing_record = None
                else:
                    version_comment = ""
                    existing_record = record.__class__.filter(
                        **{name_field: kwargs[name_field]}
                    ).one_or_none()
                if existing_record is not None:
                    logger.important(
                        f"returning existing {record.__class__.__name__} record with same"
                        f" {name_field}{version_comment}: '{kwargs[name_field]}'"
                    )
                    init_self_from_db(record, existing_record)
                    return None
        super(Record, record).__init__(**kwargs)
    elif len(args) != len(record._meta.concrete_fields):
        raise ValueError("please provide keyword arguments, not plain arguments")
    else:
        # object is loaded from DB (**kwargs could be omitted below, I believe)
        super(Record, record).__init__(*args, **kwargs)
        _store_record_old_name(record)


@classmethod  # type:ignore
@doc_args(Record.filter.__doc__)
def filter(cls, *queries, **expressions) -> QuerySet:
    """{}"""  # noqa: D415
    from lamindb._query_set import QuerySet

    _using_key = None
    if "_using_key" in expressions:
        _using_key = expressions.pop("_using_key")

    return QuerySet(model=cls, using=_using_key).filter(*queries, **expressions)


@classmethod  # type:ignore
@doc_args(Record.get.__doc__)
def get(
    cls,
    idlike: int | str | None = None,
    **expressions,
) -> Record:
    """{}"""  # noqa: D415
    from lamindb._query_set import QuerySet

    return QuerySet(model=cls).get(idlike, **expressions)


@classmethod  # type:ignore
@doc_args(Record.df.__doc__)
def df(
    cls,
    include: str | list[str] | None = None,
    join: str = "inner",
    limit: int = 100,
) -> pd.DataFrame:
    """{}"""  # noqa: D415
    from lamindb._filter import filter

    query_set = filter(cls)
    if hasattr(cls, "updated_at"):
        query_set = query_set.order_by("-updated_at")
    return query_set[:limit].df(include=include, join=join)


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
    registry = input_queryset.model
    if field is None:
        fields = [
            field.name
            for field in registry._meta.fields
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
    refined_output_queryset = output_queryset.filter(narrow_expression).annotate(
        ordering=Value(1, output_field=IntegerField())
    )
    remaining_output_queryset = output_queryset.exclude(narrow_expression).annotate(
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
    field = get_name_field(registry=queryset.model, field=field)

    return Lookup(
        records=queryset,
        values=[i.get(field) for i in queryset.values()],
        tuple_name=cls.__class__.__name__,
        prefix="ln",
    ).lookup(
        return_field=(
            get_name_field(registry=queryset.model, field=return_field)
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


def get_name_field(
    registry: type[Record] | QuerySet | Manager,
    *,
    field: str | StrField | None = None,
) -> str:
    """Get the 1st char or text field from the registry."""
    if isinstance(registry, (QuerySet, Manager)):
        registry = registry.model
    model_field_names = [i.name for i in registry._meta.fields]

    # set to default name field
    if field is None:
        if hasattr(registry, "_name_field"):
            field = registry._meta.get_field(registry._name_field)
        elif "name" in model_field_names:
            field = registry._meta.get_field("name")
        else:
            # first char or text field that doesn't contain "id"
            for i in registry._meta.fields:
                if "id" in i.name:
                    continue
                if i.get_internal_type() in {"CharField", "TextField"}:
                    field = i
                    break

        # no default name field can be found
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
    elif using_key is None or using_key == "default":
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
    instance: str | None,
) -> QuerySet:
    """{}"""  # noqa: D415
    from ._query_set import QuerySet

    if instance is None:
        return QuerySet(model=cls, using=None)
    owner, name = get_owner_name_from_identifier(instance)
    settings_file = instance_settings_file(name, owner)
    cache_filepath = ln_setup.settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
    if not settings_file.exists():
        result = connect_instance_hub(owner=owner, name=name)
        if isinstance(result, str):
            raise RuntimeError(
                f"Failed to load instance {instance}, please check your permissions!"
            )
        iresult, _ = result
        source_schema = {
            schema for schema in iresult["schema_str"].split(",") if schema != ""
        }  # type: ignore
        target_schema = ln_setup.settings.instance.schema
        if not source_schema.issubset(target_schema):
            missing_members = source_schema - target_schema
            logger.warning(
                f"source schema has additional modules: {missing_members}\nconsider mounting these schema modules to transfer all metadata"
            )
        cache_filepath.write_text(f"{iresult['lnid']}\n{iresult['schema_str']}")  # type: ignore
        settings_file = instance_settings_file(name, owner)
        db = update_db_using_local(iresult, settings_file)
    else:
        isettings = load_instance_settings(settings_file)
        db = isettings.db
        cache_filepath.write_text(f"{isettings.uid}\n{','.join(isettings.schema)}")  # type: ignore
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
    transfer_logs: dict,
):
    record = records[0] if isinstance(records, (list, QuerySet)) else records
    if hasattr(record, f"{fk}_id") and getattr(record, f"{fk}_id") is not None:
        fk_record = getattr(record, fk)
        field = REGISTRY_UNIQUE_FIELD.get(fk, "uid")
        fk_record_default = fk_record.__class__.filter(
            **{field: getattr(fk_record, field)}
        ).one_or_none()
        if fk_record_default is None:
            from copy import copy

            fk_record_default = copy(fk_record)
            transfer_to_default_db(
                fk_record_default, using_key, save=True, transfer_logs=transfer_logs
            )
        if isinstance(records, (list, QuerySet)):
            for r in records:
                setattr(r, f"{fk}", None)
                setattr(r, f"{fk}_id", fk_record_default.id)
        else:
            setattr(records, f"{fk}", None)
            setattr(records, f"{fk}_id", fk_record_default.id)


FKBULK = [
    "organism",
    "source",
    "_source_code_artifact",  # Transform
    "report",  # Run
]


def transfer_fk_to_default_db_bulk(
    records: list | QuerySet, using_key: str | None, transfer_logs: dict
):
    for fk in FKBULK:
        update_fk_to_default_db(records, fk, using_key, transfer_logs=transfer_logs)


def get_transfer_run(record) -> Run:
    from lamindb_setup import settings as setup_settings

    from lamindb.core._context import context
    from lamindb.core._data import WARNING_RUN_TRANSFORM

    slug = record._state.db
    owner, name = get_owner_name_from_identifier(slug)
    cache_filepath = ln_setup.settings.cache_dir / f"instance--{owner}--{name}--uid.txt"
    if not cache_filepath.exists():
        raise SystemExit("Need to call .using() before")
    instance_uid = cache_filepath.read_text().split("\n")[0]
    key = f"transfers/{instance_uid}"
    uid = instance_uid + "0000"
    transform = Transform.filter(uid=uid).one_or_none()
    if transform is None:
        search_names = settings.creation.search_names
        settings.creation.search_names = False
        transform = Transform(
            uid=uid, name=f"Transfer from `{slug}`", key=key, type="function"
        ).save()
        settings.creation.search_names = search_names
    # use the global run context to get the parent run id
    if context.run is not None:
        parent = context.run
    else:
        if not settings.creation.artifact_silence_missing_run_warning:
            logger.warning(WARNING_RUN_TRANSFORM)
        parent = None
    # it doesn't seem to make sense to create new runs for every transfer
    run = Run.filter(transform=transform, parent=parent).one_or_none()
    if run is None:
        run = Run(transform=transform, parent=parent).save()
        run.parent = parent  # so that it's available in memory
    return run


def transfer_to_default_db(
    record: Record,
    using_key: str | None,
    *,
    transfer_logs: dict,
    save: bool = False,
    transfer_fk: bool = True,
) -> Record | None:
    if record._state.db is None or record._state.db == "default":
        return None
    registry = record.__class__
    record_on_default = registry.objects.filter(uid=record.uid).one_or_none()
    record_str = f"{record.__class__.__name__}(uid='{record.uid}')"
    if transfer_logs["run"] is None:
        transfer_logs["run"] = get_transfer_run(record)
    if record_on_default is not None:
        transfer_logs["mapped"].append(record_str)
        return record_on_default
    else:
        transfer_logs["transferred"].append(record_str)

    if hasattr(record, "created_by_id"):
        record.created_by = None
        record.created_by_id = ln_setup.settings.user.id
    # run & transform
    run = transfer_logs["run"]
    if hasattr(record, "run_id"):
        record.run = None
        record.run_id = run.id
    # deal with denormalized transform FK on artifact and collection
    if hasattr(record, "transform_id"):
        record.transform = None
        record.transform_id = run.transform_id
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
        update_fk_to_default_db(record, fk, using_key, transfer_logs=transfer_logs)
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
        artifacts = self.ordered_artifacts.list()
    pre_existing_record = None
    # consider records that are being transferred from other databases
    transfer_logs: dict[str, list[str]] = {"mapped": [], "transferred": [], "run": None}
    if db is not None and db != "default" and using_key is None:
        if isinstance(self, IsVersioned):
            if not self.is_latest:
                raise NotImplementedError(
                    "You are attempting to transfer a record that's not the latest in its version history. This is currently not supported."
                )
        pre_existing_record = transfer_to_default_db(
            self, using_key, transfer_logs=transfer_logs
        )
    if pre_existing_record is not None:
        init_self_from_db(self, pre_existing_record)
    else:
        # save versioned record
        if isinstance(self, IsVersioned) and self._revises is not None:
            assert self._revises.is_latest  # noqa: S101
            revises = self._revises
            revises.is_latest = False
            with transaction.atomic():
                revises._revises = None  # ensure we don't start a recursion
                revises.save()
                check_name_change(self)
                super(Record, self).save(*args, **kwargs)
                _store_record_old_name(self)
            self._revises = None
        # save unversioned record
        else:
            check_name_change(self)
            super(Record, self).save(*args, **kwargs)
            _store_record_old_name(self)
    # perform transfer of many-to-many fields
    # only supported for Artifact and Collection records
    if db is not None and db != "default" and using_key is None:
        if self.__class__.__name__ == "Collection":
            if len(artifacts) > 0:
                logger.info("transfer artifacts")
                for artifact in artifacts:
                    artifact.save()
                self.artifacts.add(*artifacts)
        if hasattr(self, "labels"):
            from copy import copy

            from lnschema_core.models import FeatureManager

            # here we go back to original record on the source database
            self_on_db = copy(self)
            self_on_db._state.db = db
            self_on_db.pk = pk_on_db  # manually set the primary key
            self_on_db.features = FeatureManager(self_on_db)
            self.features._add_from(self_on_db, transfer_logs=transfer_logs)
            self.labels.add_from(self_on_db, transfer_logs=transfer_logs)
        for k, v in transfer_logs.items():
            if k != "run":
                logger.important(f"{k} records: {', '.join(v)}")
    return self


def _store_record_old_name(record: Record):
    # writes the name to the _name attribute, so we can detect renaming upon save
    if hasattr(record, "_name_field"):
        record._name = getattr(record, record._name_field)


def check_name_change(record: Record):
    """Warns if a record's name has changed."""
    if (
        not record.pk
        or not hasattr(record, "_name")
        or not hasattr(record, "_name_field")
    ):
        return

    old_name = record._name
    new_name = getattr(record, record._name_field)
    registry = record.__class__.__name__

    if old_name != new_name:
        # when a label is renamed, only raise a warning if it has a feature
        if hasattr(record, "artifacts"):
            linked_records = (
                record.artifacts.through.filter(
                    label_ref_is_name=True, **{f"{registry.lower()}_id": record.pk}
                )
                .exclude(feature_id=None)  # must have a feature
                .exclude(
                    feature_ref_is_name=None
                )  # must be linked via Curator and therefore part of a featureset
                .distinct()
            )
            artifact_ids = linked_records.list("artifact__uid")
            n = len(artifact_ids)
            s = "s" if n > 1 else ""
            if n > 0:
                logger.error(
                    f"You are trying to {colors.red('rename label')} from '{old_name}' to '{new_name}'!\n"
                    f"   → The following {n} artifact{s} {colors.red('will no longer be validated')}: {artifact_ids}\n\n"
                    f"{colors.bold('To rename this label')}, make it external:\n"
                    f"   → run `artifact.labels.make_external(label)`\n\n"
                    f"After renaming, consider re-curating the above artifact{s}:\n"
                    f'   → in each dataset, manually modify label "{old_name}" to "{new_name}"\n'
                    f"   → run `ln.Curator`\n"
                )
                raise RecordNameChangeIntegrityError

        # when a feature is renamed
        elif isinstance(record, Feature):
            # only internal features are associated with featuresets
            linked_artifacts = Artifact.filter(feature_sets__features=record).list(
                "uid"
            )
            n = len(linked_artifacts)
            s = "s" if n > 1 else ""
            if n > 0:
                logger.error(
                    f"You are trying to {colors.red('rename feature')} from '{old_name}' to '{new_name}'!\n"
                    f"   → The following {n} artifact{s} {colors.red('will no longer be validated')}: {linked_artifacts}\n\n"
                    f"{colors.bold('To rename this feature')}, make it external:\n"
                    "   → run `artifact.features.make_external(feature)`\n\n"
                    f"After renaming, consider re-curating the above artifact{s}:\n"
                    f"   → in each dataset, manually modify feature '{old_name}' to '{new_name}'\n"
                    f"   → run `ln.Curator`\n"
                )
                raise RecordNameChangeIntegrityError


def delete(self) -> None:
    """Delete the record."""
    # note that the logic below does not fire if a record is moved to the trash
    # the idea is that moving a record to the trash should move its entire version family
    # to the trash, whereas permanently deleting should default to only deleting a single record
    # of a version family
    # we can consider making it easy to permanently delete entire version families as well,
    # but that's for another time
    if isinstance(self, IsVersioned) and self.is_latest:
        new_latest = (
            self.__class__.objects.using(self._state.db)
            .filter(is_latest=False, uid__startswith=self.stem_uid)
            .order_by("-created_at")
            .first()
        )
        if new_latest is not None:
            new_latest.is_latest = True
            with transaction.atomic():
                new_latest.save()
                super(Record, self).delete()
            logger.warning(f"new latest version is {new_latest}")
            return None
    super(Record, self).delete()


METHOD_NAMES = [
    "__init__",
    "filter",
    "get",
    "df",
    "search",
    "lookup",
    "save",
    "delete",
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
