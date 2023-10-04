import builtins
from typing import Iterable, List, NamedTuple, Optional, Union

import lamindb_setup as ln_setup
import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models import Manager, QuerySet
from lamin_utils import logger
from lamin_utils._lookup import Lookup
from lamin_utils._search import search as base_search
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Registry
from lnschema_core.types import ListLike, StrField

from lamindb._utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records

IPYTHON = getattr(builtins, "__IPYTHON__", False)


class ValidationError(Exception):
    pass


def init_self_from_db(self: Registry, existing_record: Registry):
    new_args = [
        getattr(existing_record, field.attname) for field in self._meta.concrete_fields
    ]
    super(self.__class__, self).__init__(*new_args)
    self._state.adding = False  # mimic from_db
    self._state.db = "default"


def validate_required_fields(orm: Registry, kwargs):
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


def suggest_objects_with_same_name(orm: Registry, kwargs) -> Optional[str]:
    if kwargs.get("name") is None:
        return None
    else:
        results = orm.search(kwargs["name"])
        if results.shape[0] == 0:
            return None

        # subset results to those with at least 0.90 levensteihn distance
        results = results.loc[results.__ratio__ >= 90]

        # test for exact match
        if len(results) > 0:
            if results.index[0] == kwargs["name"]:
                return "object-with-same-name-exists"
            else:
                s = "" if results.shape[0] == 1 else "s"
                it = "it" if results.shape[0] == 1 else "one of them"
                msg = (
                    f"record{s} with similar name{s} exist! did you mean to load {it}?"
                )
                if IPYTHON:
                    from IPython.display import display

                    logger.warning(f"{msg}")
                    display(results)
                else:
                    logger.warning(f"{msg}\n{results}")
    return None


def __init__(orm: Registry, *args, **kwargs):
    if not args:
        validate_required_fields(orm, kwargs)
        from .dev._settings import settings

        # do not search for names if an id is passed; this is important
        # e.g. when synching ids from the notebook store to lamindb
        has_consciously_provided_id = False
        if "_has_consciously_provided_id" in kwargs:
            has_consciously_provided_id = kwargs.pop("_has_consciously_provided_id")
        if settings.upon_create_search_names and not has_consciously_provided_id:
            result = suggest_objects_with_same_name(orm, kwargs)
            if result == "object-with-same-name-exists":
                if "version" in kwargs:
                    version_comment = " and version"
                    existing_record = orm.filter(
                        name=kwargs["name"], version=kwargs["version"]
                    ).one_or_none()
                else:
                    version_comment = ""
                    existing_record = orm.filter(name=kwargs["name"]).one()
                if existing_record is not None:
                    logger.success(
                        f"loaded {orm.__class__.__name__} record with exact same"
                        f" name{version_comment}: '{kwargs['name']}'"
                    )
                    init_self_from_db(orm, existing_record)
                    return None
        super(Registry, orm).__init__(**kwargs)
    elif len(args) != len(orm._meta.concrete_fields):
        raise ValueError("please provide keyword arguments, not plain arguments")
    else:
        # object is loaded from DB (**kwargs could be omitted below, I believe)
        super(Registry, orm).__init__(*args, **kwargs)


# from_values doesn't apply for QuerySet or Manager
@classmethod  # type:ignore
@doc_args(Registry.from_values.__doc__)
def from_values(
    cls, values: ListLike, field: Optional[StrField] = None, **kwargs
) -> List["Registry"]:
    """{}"""
    from_bionty = True if cls.__module__.startswith("lnschema_bionty.") else False
    field_str = get_default_str_field(cls, field=field)
    return get_or_create_records(
        iterable=values,
        field=getattr(cls, field_str),
        from_bionty=from_bionty,
        **kwargs,
    )


# From: https://stackoverflow.com/a/37648265
def _order_queryset_by_ids(queryset: QuerySet, ids: Iterable):
    from django.db.models import Case, When

    preserved = Case(*[When(pk=pk, then=pos) for pos, pk in enumerate(ids)])
    return queryset.filter(pk__in=ids).order_by(preserved)


def _search(
    cls,
    string: str,
    *,
    field: Optional[Union[StrField, List[StrField]]] = None,
    limit: Optional[int] = 20,
    return_queryset: bool = False,
    case_sensitive: bool = False,
    synonyms_field: Optional[StrField] = "synonyms",
) -> Union["pd.DataFrame", "QuerySet"]:
    queryset = _queryset(cls)
    orm = queryset.model

    def _search_single_field(
        string: str,
        field: Optional[StrField],
        synonyms_field: Optional[StrField] = "synonyms",
    ) -> "pd.DataFrame":
        field = get_default_str_field(orm=orm, field=field)

        try:
            orm._meta.get_field(synonyms_field)
            synonyms_field_exists = True
        except FieldDoesNotExist:
            synonyms_field_exists = False

        if synonyms_field is not None and synonyms_field_exists:
            df = pd.DataFrame(queryset.values("id", field, synonyms_field))
        else:
            df = pd.DataFrame(queryset.values("id", field))

        return base_search(
            df=df,
            string=string,
            field=field,
            limit=limit,
            synonyms_field=str(synonyms_field),
            case_sensitive=case_sensitive,
        )

    # search in both key and description fields for file
    if orm._meta.model.__name__ == "File" and field is None:
        field = ["key", "description"]

    if not isinstance(field, List):
        field = [field]

    results = []
    for fd in field:
        result_field = _search_single_field(
            string=string, field=fd, synonyms_field=synonyms_field
        )
        results.append(result_field)
        # turn off synonyms search after the 1st field
        synonyms_field = None

    if len(results) > 1:
        result = (
            pd.concat([r.reset_index() for r in results], join="outer")
            .drop(columns=["index"], errors="ignore")
            .set_index("id")
        )
    else:
        result = results[0]

    # remove results that have __ratio__ 0
    if "__ratio__" in result.columns:
        result = result[result["__ratio__"] > 0].sort_values(
            "__ratio__", ascending=False
        )
        # move the __ratio__ to be the last column
        result["__ratio__"] = result.pop("__ratio__")

    if return_queryset:
        return _order_queryset_by_ids(queryset, result.reset_index()["id"])
    else:
        return result.fillna("")


@classmethod  # type: ignore
@doc_args(Registry.search.__doc__)
def search(
    cls,
    string: str,
    *,
    field: Optional[StrField] = None,
    limit: Optional[int] = 20,
    return_queryset: bool = False,
    case_sensitive: bool = False,
    synonyms_field: Optional[StrField] = "synonyms",
) -> Union["pd.DataFrame", "QuerySet"]:
    """{}"""
    return _search(
        cls=cls,
        string=string,
        field=field,
        return_queryset=return_queryset,
        limit=limit,
        case_sensitive=case_sensitive,
        synonyms_field=synonyms_field,
    )


def _lookup(
    cls, field: Optional[StrField] = None, return_field: Optional[StrField] = None
) -> NamedTuple:
    """{}"""
    queryset = _queryset(cls)
    field = get_default_str_field(orm=queryset.model, field=field)

    return Lookup(
        records=queryset,
        values=[i.get(field) for i in queryset.values()],
        tuple_name=cls.__class__.__name__,
        prefix="ln",
    ).lookup(
        return_field=get_default_str_field(orm=queryset.model, field=return_field)
        if return_field is not None
        else None
    )


@classmethod  # type: ignore
@doc_args(Registry.lookup.__doc__)
def lookup(
    cls, field: Optional[StrField] = None, return_field: Optional[StrField] = None
) -> NamedTuple:
    """{}"""
    return _lookup(cls=cls, field=field, return_field=return_field)


def get_default_str_field(
    orm: Union[Registry, QuerySet, Manager],
    *,
    field: Optional[Union[str, StrField]] = None,
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
                "please pass a Registry string field, e.g., `CellType.name`!"
            )
        else:
            field = field.name  # type:ignore
    if not isinstance(field, str):
        try:
            field = field.field.name
        except AttributeError:
            raise TypeError(
                "please pass a Registry string field, e.g., `CellType.name`!"
            )

    return field


def _queryset(cls: Union[Registry, QuerySet, Manager]) -> QuerySet:
    queryset = cls.all() if isinstance(cls, QuerySet) else cls.objects.all()
    return queryset


def transfer_to_default_db(record: Registry, save: bool = False):
    db = record._state.db
    if db is not None and db != "default":
        logger.info(f"saving from instance {db} to default instance: {record}")
        from lamindb.dev._data import WARNING_RUN_TRANSFORM
        from lamindb.dev._run_context import run_context

        logger.hint("saving to default instance")
        if (
            hasattr(record, "created_by_id")
            and record.created_by_id != ln_setup.settings.user.id
        ):
            logger.info(f"updating created_by_id with {ln_setup.settings.user.id}")
            record.created_by_id = ln_setup.settings.user.id
        if hasattr(record, "run_id"):
            if run_context.run is not None:
                logger.info("updating run & transform to current run & transform")
                record.run_id = run_context.run.id
            else:
                logger.warning(WARNING_RUN_TRANSFORM)
                record.run_id = None
        if hasattr(record, "transform_id"):
            if run_context.transform is not None:
                record.transform_id = run_context.transform.id
            else:
                record.transform_id = None
        if hasattr(record, "storage_id") and record.storage_id is not None:
            record.storage.save()
        record._state.db = "default"
        if save:
            record.save()


# docstring handled through attach_func_to_class_method
def save(self, *args, **kwargs) -> None:
    db = self._state.db
    transfer_to_default_db(self)
    super(Registry, self).save(*args, **kwargs)
    if db is not None and db != "default":
        if hasattr(self, "labels"):
            logger.info("transfer labels")
            from copy import copy

            self_on_db = copy(self)
            self_on_db._state.db = db
            self.labels.add_from(self_on_db)


METHOD_NAMES = [
    "__init__",
    "search",
    "lookup",
    "save",
    "from_values",
]

if _TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(Registry, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Registry, globals())


@classmethod  # type: ignore
def __get_schema_name__(cls) -> str:
    schema_module_name = cls.__module__.split(".")[0]
    schema_name = schema_module_name.replace("lnschema_", "")
    return schema_name


@classmethod  # type: ignore
def __get_name_with_schema__(cls) -> str:
    schema_name = cls.__get_schema_name__()
    return f"{schema_name}.{cls.__name__}"


setattr(Registry, "__get_schema_name__", __get_schema_name__)
setattr(Registry, "__get_name_with_schema__", __get_name_with_schema__)
