import builtins
from datetime import datetime
from typing import Any, Dict, Iterable, List, Literal, NamedTuple, Optional, Set, Union

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db import models
from django.db.models import CharField, TextField
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_logger import logger
from lamin_logger._lookup import Lookup
from lamin_logger._search import search as base_search
from lnschema_core import ORM
from lnschema_core.types import ListLike, StrField

from . import _TESTING
from ._from_values import get_or_create_records
from .dev._settings import settings

_is_ipython = getattr(builtins, "__IPYTHON__", False)


class ValidationError(Exception):
    pass


def validate_required_fields(orm: ORM, kwargs):
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


def suggest_objects_with_same_name(orm: ORM, kwargs) -> Optional[str]:
    if kwargs.get("name") is None:
        return None
    else:
        results = orm.search(kwargs["name"])
        if results.shape[0] == 0:
            return None

        # subset results to those with at least 0.5 levensteihn distance
        results = results.loc[results.__ratio__ >= 90]

        # test for exact match
        if len(results) > 0:
            if results.index[0] == kwargs["name"]:
                logger.warning("Object with exact same name exists, returning it")
                return "object-with-same-name-exists"
            else:
                msg = "Entries with similar names exist:"
                if _is_ipython:
                    from IPython.display import display

                    logger.warning(f"{msg}")
                    display(results)
                else:
                    logger.warning(f"{msg}\n{results.name}")
    return None


def __init__(orm: ORM, *args, **kwargs):
    if not args:
        validate_required_fields(orm, kwargs)
        if settings.upon_create_search_names:
            result = suggest_objects_with_same_name(orm, kwargs)
            if result == "object-with-same-name-exists":
                existing_object = orm.select(name=kwargs["name"])[0]
                new_args = [
                    getattr(existing_object, field.attname)
                    for field in orm._meta.concrete_fields
                ]
                super(ORM, orm).__init__(*new_args)
                orm._state.adding = False  # mimic from_db
                orm._state.db = "default"
                return None
        super(ORM, orm).__init__(**kwargs)
    elif len(args) != len(orm._meta.concrete_fields):
        raise ValueError("Please provide keyword arguments, not plain arguments")
    else:
        # object is loaded from DB (**kwargs could be omitted below, I believe)
        super(ORM, orm).__init__(*args, **kwargs)


@classmethod  # type:ignore
def from_values(cls, identifiers: ListLike, field: StrField, **kwargs):
    if isinstance(field, str):
        field = getattr(cls, field)
    if not isinstance(field, Field):  # field is DeferredAttribute
        raise TypeError(
            "field must be a string or an ORM field, e.g., `CellType.name`!"
        )
    from_bionty = True if cls.__module__.startswith("lnschema_bionty.") else False
    return get_or_create_records(
        iterable=identifiers, field=field, from_bionty=from_bionty, **kwargs
    )


@classmethod  # type: ignore
def search(
    cls,
    string: str,
    *,
    field: Optional[StrField] = None,
    top_hit: bool = False,
    case_sensitive: bool = True,
    synonyms_field: Optional[Union[str, TextField, CharField]] = "synonyms",
    synonyms_sep: str = "|",
) -> Union["pd.DataFrame", "ORM"]:
    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.objects.all()
    df = pd.DataFrame.from_records(records.values())

    result = base_search(
        df=df,
        string=string,
        field=field,
        synonyms_field=str(synonyms_field),
        case_sensitive=case_sensitive,
        return_ranked_results=not top_hit,
        synonyms_sep=synonyms_sep,
        tuple_name=cls.__name__,
    )

    if not top_hit or result is None:
        return result
    else:
        if isinstance(result, list):
            return [records.get(id=r.id) for r in result]
        else:
            return records.get(id=result.id)


@classmethod  # type: ignore
def lookup(cls, field: Optional[StrField] = None) -> NamedTuple:
    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.objects.all()

    return Lookup(
        records=records,
        values=[i.get(field) for i in records.values()],
        tuple_name=cls.__name__,
        prefix="ln",
    ).lookup()


@classmethod  # type: ignore
def inspect(
    cls,
    identifiers: Iterable,
    field: StrField,
    *,
    case_sensitive: bool = False,
    inspect_synonyms: bool = True,
    return_df: bool = False,
    logging: bool = True,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    from lamin_logger._inspect import inspect

    if not isinstance(field, str):
        field = field.field.name

    return inspect(
        df=_filter_df_based_on_species(orm=cls, species=kwargs.get("species")),
        identifiers=identifiers,
        field=str(field),
        case_sensitive=case_sensitive,
        inspect_synonyms=inspect_synonyms,
        return_df=return_df,
        logging=logging,
    )


@classmethod  # type: ignore
def map_synonyms(
    cls,
    synonyms: Iterable,
    *,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    synonyms_sep: str = "|",
    field: Optional[str] = None,
    **kwargs,
) -> Union[List[str], Dict[str, str]]:
    from lamin_logger._map_synonyms import map_synonyms

    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    try:
        cls._meta.get_field(synonyms_field)
        df = _filter_df_based_on_species(orm=cls, species=kwargs.get("species"))
    except FieldDoesNotExist:
        df = pd.DataFrame()
    return map_synonyms(
        df=df,
        identifiers=synonyms,
        field=field,
        return_mapper=return_mapper,
        case_sensitive=case_sensitive,
        keep=keep,
        synonyms_field=synonyms_field,
        sep=synonyms_sep,
    )


def _filter_df_based_on_species(orm: ORM, species: Union[str, ORM, None] = None):
    import pandas as pd

    records = orm.objects.all()
    try:
        # if the orm has a species field, it's required
        records.model._meta.get_field("species")
        if species is None:
            raise AssertionError(
                f"{orm.__name__} table requires to specify a species name via"
                " `species=`!"
            )
        elif isinstance(species, ORM):
            species_name = species.name
        else:
            species_name = species
        records = records.filter(species__name=species_name)
    except FieldDoesNotExist:
        pass

    return pd.DataFrame.from_records(records.values())


def get_default_str_field(orm: ORM) -> str:
    """Get the 1st char or text field from the orm."""
    model_field_names = [i.name for i in orm._meta.fields]

    # set default field
    if "name" in model_field_names:
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
        raise ValueError("Please specify a field to search against!")

    return field.name


def _add_or_remove_synonyms(
    synonym: Union[str, Iterable],
    record: ORM,
    action: Literal["add", "remove"],
    force: bool = False,
):
    """Add or remove synonyms."""

    def check_synonyms_in_all_records(synonyms: Set[str], record: ORM):
        """Errors if input synonyms are already associated with records in the DB."""
        import pandas as pd
        from IPython.display import display

        syns_all = (
            record.__class__.objects.exclude(synonyms="").exclude(synonyms=None).all()
        )
        if len(syns_all) == 0:
            return
        df = pd.DataFrame(syns_all.values())
        df["synonyms"] = df["synonyms"].str.split("|")
        df = df.explode("synonyms")
        matches_df = df[(df["synonyms"].isin(synonyms)) & (df["id"] != record.id)]
        if matches_df.shape[0] > 0:
            records_df = pd.DataFrame(syns_all.filter(id__in=matches_df["id"]).values())
            logger.error(
                f"Input synonyms {matches_df['synonyms'].unique()} already associated"
                " with the following records:\n(Pass `force=True` to ignore this error)"
            )
            display(records_df)
            raise SystemExit(AssertionError)

    # passed synonyms
    if isinstance(synonym, str):
        syn_new_set = set([synonym])
    else:
        syn_new_set = set(synonym)
    # nothing happens when passing an empty string or list
    if len(syn_new_set) == 0:
        return
    # because we use | as the separator
    if any(["|" in i for i in syn_new_set]):
        raise AssertionError("A synonym can't contain '|'!")

    # existing synonyms
    syns_exist = record.synonyms
    if syns_exist is None or len(syns_exist) == 0:
        syns_exist_set = set()
    else:
        syns_exist_set = set(syns_exist.split("|"))

    if action == "add":
        if not force:
            check_synonyms_in_all_records(syn_new_set, record)
        syns_exist_set.update(syn_new_set)
    elif action == "remove":
        syns_exist_set = syns_exist_set.difference(syn_new_set)

    if len(syns_exist_set) == 0:
        syns_str = None
    else:
        syns_str = "|".join(syns_exist_set)

    record.synonyms = syns_str

    # if the record already exists in the DB, save it
    if not record._state.adding:
        record.save()


def _check_synonyms_field_exist(record: ORM):
    try:
        record.__getattribute__("synonyms")
    except AttributeError:
        raise NotImplementedError(
            f"No synonyms field found in table {record.__class__.__name__}!"
        )


def add_synonym(self, synonym: Union[str, ListLike], force: bool = False):
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(synonym=synonym, record=self, force=force, action="add")


def remove_synonym(self, synonym: Union[str, ListLike]):
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(synonym=synonym, record=self, action="remove")


def format_datetime(dt: Union[datetime, Any]) -> str:
    if not isinstance(dt, datetime):
        return dt
    else:
        return dt.strftime("%Y-%m-%d %H:%M:%S")


def __repr__(self: ORM) -> str:
    field_names = [
        field.name
        for field in self._meta.fields
        if not isinstance(field, (models.ForeignKey, models.DateTimeField))
    ]
    # skip created_at
    field_names += [
        field.name
        for field in self._meta.fields
        if isinstance(field, models.DateTimeField) and field.name != "created_at"
    ]
    field_names += [
        f"{field.name}_id"
        for field in self._meta.fields
        if isinstance(field, models.ForeignKey)
    ]
    fields_str = {
        k: format_datetime(getattr(self, k)) for k in field_names if hasattr(self, k)
    }
    fields_joined_str = ", ".join([f"{k}={fields_str[k]}" for k in fields_str])
    return f"{self.__class__.__name__}({fields_joined_str})"


# this captures the original signatures for testing purposes
# it's used in the unit tests
if _TESTING:
    from inspect import signature

    SIG_ORM_SEARCH = signature(ORM.search)
    SIG_ORM_LOOKUP = signature(ORM.lookup)
    SIG_ORM_INSPECT = signature(ORM.inspect)
    SIG_ORM_FROM_VALUES = signature(ORM.from_values)
    SIG_ORM_MAP_SYNONYM = signature(ORM.map_synonyms)
    SIG_ORM_ADD_SYNONYM = signature(ORM.add_synonym)
    SIG_ORM_REMOVE_SYNONYM = signature(ORM.remove_synonym)


ORM.__init__ = __init__
ORM.__str__ = __repr__
ORM.search = search
ORM.lookup = lookup
ORM.map_synonyms = map_synonyms
ORM.inspect = inspect
ORM.add_synonym = add_synonym
ORM.remove_synonym = remove_synonym
ORM.from_values = from_values
