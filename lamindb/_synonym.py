from typing import Dict, Iterable, List, Literal, Optional, Set, Union

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models import QuerySet
from lamin_utils import logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Registry, SynonymsAware
from lnschema_core.types import ListLike

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._registry import get_default_str_field
from ._validate import _filter_query_based_on_species


@classmethod  # type: ignore
@doc_args(SynonymsAware.map_synonyms.__doc__)
def map_synonyms(
    cls,
    synonyms: Iterable,
    *,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    field: Optional[str] = None,
    **kwargs,
) -> Union[List[str], Dict[str, str]]:
    """{}"""
    return _map_synonyms(
        cls=cls,
        synonyms=synonyms,
        return_mapper=return_mapper,
        case_sensitive=case_sensitive,
        keep=keep,
        synonyms_field=synonyms_field,
        field=field,
        **kwargs,
    )


def set_abbr(self, value: str):
    try:
        self.add_synonym(value, save=False)
    except NotImplementedError:
        pass
    self.abbr = value
    if not self._state.adding:
        self.save()


def add_synonym(
    self,
    synonym: Union[str, ListLike],
    force: bool = False,
    save: Optional[bool] = None,
):
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(
        synonym=synonym, record=self, force=force, action="add", save=save
    )


def remove_synonym(self, synonym: Union[str, ListLike]):
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(synonym=synonym, record=self, action="remove")


def _add_or_remove_synonyms(
    synonym: Union[str, Iterable],
    record: Registry,
    action: Literal["add", "remove"],
    force: bool = False,
    save: Optional[bool] = None,
):
    """Add or remove synonyms."""

    def check_synonyms_in_all_records(synonyms: Set[str], record: Registry):
        """Errors if input synonym is associated with other records in the DB."""
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
                f"input synonyms {matches_df['synonyms'].unique()} already associated"
                " with the following records:\n"
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
        raise AssertionError("a synonym can't contain '|'!")

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

    if save is None:
        # if record is already in DB, save the changes to DB
        save = not record._state.adding
    if save:
        record.save()


def _check_synonyms_field_exist(record: Registry):
    try:
        record.__getattribute__("synonyms")
    except AttributeError:
        raise NotImplementedError(
            f"No synonyms field found in table {record.__class__.__name__}!"
        )


def _map_synonyms(
    cls,
    synonyms: Iterable,
    *,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    field: Optional[str] = None,
    **kwargs,
) -> Union[List[str], Dict[str, str]]:
    """{}"""
    from lamin_utils._map_synonyms import map_synonyms

    if isinstance(synonyms, str):
        synonyms = [synonyms]
    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    cls = cls.model if isinstance(cls, QuerySet) else cls

    try:
        cls._meta.get_field(synonyms_field)
        df = _filter_query_based_on_species(orm=cls, species=kwargs.get("species"))
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
    )


METHOD_NAMES = ["map_synonyms", "add_synonym", "remove_synonym", "set_abbr"]

if _TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(SynonymsAware, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, SynonymsAware, globals())
