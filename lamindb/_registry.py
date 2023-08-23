import builtins
from typing import Iterable, List, NamedTuple, Optional, Union

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models import Manager, QuerySet
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_utils import logger
from lamin_utils._lookup import Lookup
from lamin_utils._search import search as base_search
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Registry
from lnschema_core.types import ListLike, StrField

from lamindb.dev.utils import attach_func_to_class_method

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

        # subset results to those with at least 0.85 levensteihn distance
        results = results.loc[results.__ratio__ >= 85]

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

        if settings.upon_create_search_names:
            result = suggest_objects_with_same_name(orm, kwargs)
            if result == "object-with-same-name-exists":
                if "version" in kwargs:
                    version_comment = " and version "
                    existing_record = orm.filter(
                        name=kwargs["name"], version=kwargs["version"]
                    ).one_or_none()
                else:
                    version_comment = " "
                    existing_record = orm.filter(name=kwargs["name"]).one()
                if existing_record is not None:
                    logger.success(
                        f"loaded record with exact same name{version_comment}"
                    )
                    init_self_from_db(orm, existing_record)
                    return None
        super(Registry, orm).__init__(**kwargs)
    elif len(args) != len(orm._meta.concrete_fields):
        raise ValueError("please provide keyword arguments, not plain arguments")
    else:
        # object is loaded from DB (**kwargs could be omitted below, I believe)
        super(Registry, orm).__init__(*args, **kwargs)


@classmethod  # type:ignore
@doc_args(Registry.from_values.__doc__)
def from_values(cls, values: ListLike, field: StrField, **kwargs) -> List["Registry"]:
    """{}"""
    if isinstance(field, str):
        field = getattr(cls, field)
    if not isinstance(field, Field):  # field is DeferredAttribute
        raise TypeError(
            "field must be a string or an Registry field, e.g., `CellType.name`!"
        )
    from_bionty = True if cls.__module__.startswith("lnschema_bionty.") else False
    return get_or_create_records(
        iterable=values, field=field, from_bionty=from_bionty, **kwargs
    )


# From: https://stackoverflow.com/a/37648265
def _order_query_set_by_ids(queryset: QuerySet, ids: Iterable):
    from django.db.models import Case, When

    preserved = Case(*[When(pk=pk, then=pos) for pos, pk in enumerate(ids)])
    return queryset.filter(pk__in=ids).order_by(preserved)


def _search(
    cls,
    string: str,
    *,
    field: Optional[Union[StrField, List[StrField]]] = None,
    return_queryset: bool = False,
    limit: Optional[int] = None,
    case_sensitive: bool = False,
    synonyms_field: Optional[StrField] = "synonyms",
) -> Union["pd.DataFrame", "QuerySet"]:
    query_set = cls.all() if isinstance(cls, QuerySet) else cls.objects.all()
    orm = cls.model if isinstance(cls, QuerySet) else cls

    def _search_single_field(
        string: str,
        field: Optional[StrField],
        synonyms_field: Optional[StrField] = "synonyms",
    ) -> "pd.DataFrame":
        if field is None:
            field = get_default_str_field(cls)
        if not isinstance(field, str):
            field = field.field.name

        try:
            orm._meta.get_field(synonyms_field)
            synonyms_field_exists = True
        except FieldDoesNotExist:
            synonyms_field_exists = False

        if synonyms_field is not None and synonyms_field_exists:
            df = pd.DataFrame(query_set.values("id", field, synonyms_field))
        else:
            df = pd.DataFrame(query_set.values("id", field))

        return base_search(
            df=df,
            string=string,
            field=field,
            limit=limit,
            synonyms_field=str(synonyms_field),
            case_sensitive=case_sensitive,
        )

    # search in both key and description fields for file
    if orm.__name__ == "File" and field is None:
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
        return _order_query_set_by_ids(query_set, result.reset_index()["id"])
    else:
        return result.fillna("")


@classmethod  # type: ignore
@doc_args(Registry.search.__doc__)
def search(
    cls,
    string: str,
    *,
    field: Optional[StrField] = None,
    return_queryset: bool = False,
    limit: Optional[int] = None,
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


def _lookup(cls, field: Optional[StrField] = None) -> NamedTuple:
    """{}"""
    if field is None:
        if cls._meta.model.__name__ == "User":
            field = cls._meta.get_field("handle").name
        else:
            field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.all() if isinstance(cls, QuerySet) else cls.objects.all()
    cls = cls.model if isinstance(cls, QuerySet) else cls

    return Lookup(
        records=records,
        values=[i.get(field) for i in records.values()],
        tuple_name=cls.__name__,
        prefix="ln",
    ).lookup()


@classmethod  # type: ignore
@doc_args(Registry.lookup.__doc__)
def lookup(cls, field: Optional[StrField] = None) -> NamedTuple:
    """{}"""
    return _lookup(cls=cls, field=field)


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
            raise ValueError("Please specify a field!")
        else:
            field = field.name  # type:ignore
    if not isinstance(field, str):
        field = field.field.name

    return field


METHOD_NAMES = [
    "__init__",
    "search",
    "lookup",
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


def select_backward(cls, **expressions):
    logger.warning("select() is deprecated! please use: Registry.filter()")
    return cls.filter(**expressions)


@classmethod  # type: ignore
def select(cls, **expressions):
    return select_backward(cls, **expressions)


setattr(Registry, "__get_schema_name__", __get_schema_name__)
setattr(Registry, "__get_name_with_schema__", __get_name_with_schema__)
setattr(Registry, "select", select)  # backward compat
