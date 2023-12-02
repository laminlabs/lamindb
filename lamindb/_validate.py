from typing import Dict, Iterable, List, Literal, Optional, Set, Union

import numpy as np
import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models import QuerySet
from lamin_utils import colors, logger
from lamin_utils._inspect import InspectResult
from lamindb_setup.dev._docs import doc_args
from lnschema_core import CanValidate, Registry
from lnschema_core.types import ListLike, StrField

from lamindb._utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import _has_organism_field, _print_values
from ._registry import _queryset, get_default_str_field


@classmethod  # type: ignore
@doc_args(CanValidate.inspect.__doc__)
def inspect(
    cls,
    values: ListLike,
    field: Optional[Union[str, StrField]] = None,
    *,
    mute: bool = False,
    **kwargs,
) -> InspectResult:
    """{}"""
    return _inspect(
        cls=cls,
        values=values,
        field=field,
        mute=mute,
        **kwargs,
    )


@classmethod  # type: ignore
@doc_args(CanValidate.validate.__doc__)
def validate(
    cls,
    values: ListLike,
    field: Optional[Union[str, StrField]] = None,
    *,
    mute: bool = False,
    **kwargs,
) -> np.ndarray:
    """{}"""
    return _validate(cls=cls, values=values, field=field, mute=mute, **kwargs)


def _inspect(
    cls,
    values: ListLike,
    field: Optional[Union[str, StrField]] = None,
    *,
    mute: bool = False,
    **kwargs,
) -> Union["pd.DataFrame", Dict[str, List[str]]]:
    """{}"""
    from lamin_utils._inspect import inspect

    if isinstance(values, str):
        values = [values]

    field = get_default_str_field(cls, field=field)
    queryset = _queryset(cls)
    orm = queryset.model
    model_name = orm._meta.model.__name__

    # inspect in the DB
    result_db = inspect(
        df=_filter_query_based_on_organism(
            queryset=queryset, organism=kwargs.get("organism")
        ),
        identifiers=values,
        field=field,
        mute=mute,
        **kwargs,
    )
    nonval = set(result_db.non_validated).difference(result_db.synonyms_mapper.keys())

    if len(nonval) > 0 and orm.__get_schema_name__() == "bionty":
        try:
            bionty_result = orm.bionty(organism=kwargs.get("organism")).inspect(
                values=nonval, field=field, mute=True, **kwargs
            )
            bionty_validated = bionty_result.validated
            bionty_mapper = bionty_result.synonyms_mapper
            hint = False
            if len(bionty_validated) > 0 and not mute:
                print_values = _print_values(bionty_validated)
                s = "" if len(bionty_validated) == 1 else "s"
                labels = colors.yellow(f"{len(bionty_validated)} {model_name} term{s}")
                logger.print(
                    f"   detected {labels} in Bionty for"
                    f" {colors.italic(field)}: {colors.yellow(print_values)}"
                )
                hint = True

            if len(bionty_mapper) > 0 and not mute:
                print_values = _print_values(list(bionty_mapper.keys()))
                s = "" if len(bionty_mapper) == 1 else "s"
                labels = colors.yellow(f"{len(bionty_mapper)} {model_name} term{s}")
                logger.print(
                    f"   detected {labels} in Bionty as {colors.italic(f'synonym{s}')}:"
                    f" {colors.yellow(print_values)}"
                )
                hint = True

            if hint:
                logger.print(
                    f"→  add records from Bionty to your {model_name} registry via"
                    f" {colors.italic('.from_values()')}"
                )

            nonval = bionty_result.non_validated
        # no bionty source is found
        except ValueError:
            logger.warning("no Bionty source found, skipping Bionty validation")

    if len(nonval) > 0 and not mute:
        print_values = _print_values(list(nonval))
        s = "" if len(nonval) == 1 else "s"
        labels = colors.red(f"{len(nonval)} term{s}")
        logger.print(f"   couldn't validate {labels}: {colors.red(print_values)}")
        logger.print(
            f"→  if you are sure, create new record{s} via"
            f" {colors.italic(f'ln.{orm.__name__}()')} and save to your registry"
        )

    return result_db


def _validate(
    cls,
    values: ListLike,
    field: Optional[Union[str, StrField]] = None,
    *,
    mute: bool = False,
    **kwargs,
) -> np.ndarray:
    """{}"""
    from lamin_utils._inspect import validate

    return_str = True if isinstance(values, str) else False
    if isinstance(values, str):
        values = [values]

    field = get_default_str_field(cls, field=field)

    queryset = _queryset(cls)
    field_values = pd.Series(
        _filter_query_based_on_organism(
            queryset=queryset,
            organism=kwargs.get("organism"),
            values_list_field=field,
        ),
        dtype="object",
    )

    result = validate(
        identifiers=values,
        field_values=field_values,
        case_sensitive=True,
        mute=mute,
        field=field,
        **kwargs,
    )
    if return_str and len(result) == 1:
        return result[0]
    else:
        return result


@classmethod  # type: ignore
@doc_args(CanValidate.standardize.__doc__)
def standardize(
    cls,
    values: Iterable,
    field: Optional[Union[str, StrField]] = None,
    *,
    return_field: str = None,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    mute: bool = False,
    bionty_aware: bool = True,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    **kwargs,
) -> Union[List[str], Dict[str, str]]:
    """{}"""
    return _standardize(
        cls=cls,
        values=values,
        field=field,
        return_field=return_field,
        return_mapper=return_mapper,
        case_sensitive=case_sensitive,
        mute=mute,
        bionty_aware=bionty_aware,
        keep=keep,
        synonyms_field=synonyms_field,
        **kwargs,
    )


def set_abbr(self, value: str):
    self.abbr = value

    if hasattr(self, "name") and value == self.name:
        pass
    else:
        try:
            self.add_synonym(value, save=False)
        except Exception:  # pragma: no cover
            pass

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


def _standardize(
    cls,
    values: Iterable,
    field: Optional[Union[str, StrField]] = None,
    *,
    return_field: str = None,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    mute: bool = False,
    bionty_aware: bool = True,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    **kwargs,
) -> Union[List[str], Dict[str, str]]:
    """{}"""
    from lamin_utils._standardize import standardize as map_synonyms

    return_str = True if isinstance(values, str) else False
    if isinstance(values, str):
        values = [values]

    field = get_default_str_field(cls, field=field)
    return_field = get_default_str_field(
        cls, field=field if return_field is None else return_field
    )
    queryset = _queryset(cls)
    orm = queryset.model

    organism = kwargs.get("organism")
    if _has_organism_field(orm):
        # here, we can safely import lnschema_bionty
        from lnschema_bionty._bionty import create_or_get_organism_record

        organism_record = create_or_get_organism_record(organism=organism, orm=orm)
        organism = (
            organism_record.name if organism_record is not None else organism_record
        )

    try:
        orm._meta.get_field(synonyms_field)
        df = _filter_query_based_on_organism(queryset=queryset, organism=organism)
    except FieldDoesNotExist:
        df = pd.DataFrame()

    _kwargs = dict(
        field=field,
        return_field=return_field,
        case_sensitive=case_sensitive,
        keep=keep,
        synonyms_field=synonyms_field,
    )
    # standardized names from the DB
    std_names_db = map_synonyms(
        df=df,
        identifiers=values,
        return_mapper=return_mapper,
        mute=mute,
        **_kwargs,
    )

    def _return(result: List, mapper: Dict):
        if return_mapper:
            return mapper
        else:
            if return_str and len(result) == 1:
                return result[0]
            return result

    # map synonyms in Bionty
    if orm.__get_schema_name__() == "bionty" and bionty_aware:
        mapper = {}
        if return_mapper:
            mapper = std_names_db
            std_names_db = map_synonyms(
                df=df, identifiers=values, return_mapper=False, mute=True, **_kwargs
            )

        val_res = orm.validate(std_names_db, field=field, mute=True, organism=organism)
        if all(val_res):
            return _return(result=std_names_db, mapper=mapper)

        nonval = np.array(std_names_db)[~val_res]
        std_names_bt_mapper = orm.bionty(organism=organism).standardize(
            nonval, return_mapper=True, mute=True, **_kwargs
        )

        if len(std_names_bt_mapper) > 0 and not mute:
            s = "" if len(std_names_bt_mapper) == 1 else "s"
            field_print = "synonym" if field == return_field else field
            warn_msg = (
                f"found {len(std_names_bt_mapper)} {field_print}{s} in Bionty:"
                f" {list(std_names_bt_mapper.keys())}"
            )
            warn_msg += (
                f"\n   please add corresponding {orm._meta.model.__name__} records via"
                f" `.from_values({list(std_names_bt_mapper.values())})`"
            )
            logger.warning(warn_msg)

        mapper.update(std_names_bt_mapper)
        result = pd.Series(std_names_db).replace(std_names_bt_mapper).tolist()
        return _return(result=result, mapper=mapper)

    else:
        return _return(result=std_names_db, mapper=std_names_db)


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
    # nothing happens when passing an empty string or list
    if isinstance(synonym, str):
        if len(synonym) == 0:
            return
        syn_new_set = set([synonym])
    else:
        if synonym == [""]:
            return
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


def _filter_query_based_on_organism(
    queryset: QuerySet,
    organism: Optional[Union[str, Registry]] = None,
    values_list_field: Optional[str] = None,
):
    """Filter a queryset based on organism."""
    import pandas as pd

    orm = queryset.model

    if _has_organism_field(orm):
        # here, we can safely import lnschema_bionty
        from lnschema_bionty._bionty import create_or_get_organism_record

        organism_record = create_or_get_organism_record(organism=organism, orm=orm)
        if organism_record is not None:
            queryset = queryset.filter(organism__name=organism_record.name)

    if values_list_field is None:
        return pd.DataFrame.from_records(queryset.values())
    else:
        return queryset.values_list(values_list_field, flat=True)


METHOD_NAMES = [
    "validate",
    "inspect",
    "standardize",
    "add_synonym",
    "remove_synonym",
    "set_abbr",
]

if _TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(CanValidate, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, CanValidate, globals())
