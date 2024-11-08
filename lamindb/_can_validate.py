from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import lamindb_setup as ln_setup
import numpy as np
import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from lamin_utils import colors, logger
from lamindb_setup.core._docs import doc_args
from lnschema_core import CanValidate, Record

from ._from_values import _has_organism_field, _print_values, get_or_create_records
from ._record import _queryset, get_name_field
from ._utils import attach_func_to_class_method
from .core.exceptions import ValidationError

if TYPE_CHECKING:
    from django.db.models import QuerySet
    from lamin_utils._inspect import InspectResult
    from lnschema_core.types import ListLike, StrField


# from_values doesn't apply for QuerySet or Manager
@classmethod  # type:ignore
@doc_args(CanValidate.from_values.__doc__)
def from_values(
    cls,
    values: ListLike,
    field: StrField | None = None,
    create: bool = False,
    organism: Record | str | None = None,
    source: Record | None = None,
    mute: bool = False,
) -> list[Record]:
    """{}"""  # noqa: D415
    from_source = True if cls.__module__.startswith("bionty.") else False

    field_str = get_name_field(cls, field=field)
    return get_or_create_records(
        iterable=values,
        field=getattr(cls, field_str),
        create=create,
        from_source=from_source,
        organism=organism,
        source=source,
        mute=mute,
    )


@classmethod  # type: ignore
@doc_args(CanValidate.inspect.__doc__)
def inspect(
    cls,
    values: ListLike,
    field: str | StrField | None = None,
    *,
    mute: bool = False,
    organism: str | Record | None = None,
    source: Record | None = None,
) -> InspectResult:
    """{}"""  # noqa: D415
    return _inspect(
        cls=cls,
        values=values,
        field=field,
        mute=mute,
        organism=organism,
        source=source,
    )


@classmethod  # type: ignore
@doc_args(CanValidate.validate.__doc__)
def validate(
    cls,
    values: ListLike,
    field: str | StrField | None = None,
    *,
    mute: bool = False,
    organism: str | Record | None = None,
    source: Record | None = None,
) -> np.ndarray:
    """{}"""  # noqa: D415
    return _validate(
        cls=cls, values=values, field=field, mute=mute, organism=organism, source=source
    )


def _check_source_db(source: Record, using_key: str | None):
    """Check if the source is from the DB."""
    if using_key is not None and using_key != "default":
        if source._state.db != using_key:
            raise ValueError(
                f"source must be a bionty.Source record from instance '{using_key}'!"
            )


def _check_organism_db(organism: Record, using_key: str | None):
    """Check if the organism is from the DB."""
    if isinstance(organism, Record):
        if using_key is not None and using_key != "default":
            if organism._state.db != using_key:
                raise ValueError(
                    f"organism must be a bionty.Organism record from instance '{using_key}'!"
                )


def _concat_lists(values: ListLike) -> list[str]:
    """Concatenate a list of lists of strings into a single list."""
    if len(values) > 0 and isinstance(values, (list, pd.Series)):
        try:
            if isinstance(values[0], list):
                if isinstance(values, pd.Series):
                    values = values.tolist()
                values = sum([v for v in values if isinstance(v, list)], [])
        except KeyError:
            pass
    return values


def _inspect(
    cls,
    values: ListLike,
    field: str | StrField | None = None,
    *,
    mute: bool = False,
    using_key: str | None = None,
    organism: str | Record | None = None,
    source: Record | None = None,
) -> pd.DataFrame | dict[str, list[str]]:
    """{}"""  # noqa: D415
    from lamin_utils._inspect import inspect

    if isinstance(values, str):
        values = [values]
    values = _concat_lists(values)

    field = get_name_field(cls, field=field)
    queryset = _queryset(cls, using_key)
    using_key = queryset.db
    if isinstance(source, Record):
        _check_source_db(source, using_key)
        queryset = queryset.filter(source=source).all()
    _check_organism_db(organism, using_key)
    registry = queryset.model
    model_name = registry._meta.model.__name__

    # inspect in the DB
    result_db = inspect(
        df=_filter_query_based_on_organism(
            queryset=queryset, field=field, organism=organism
        ),
        identifiers=values,
        field=field,
        mute=mute,
    )
    nonval = set(result_db.non_validated).difference(result_db.synonyms_mapper.keys())

    if len(nonval) > 0 and registry.__get_schema_name__() == "bionty":
        try:
            bionty_result = registry.public(organism=organism, source=source).inspect(
                values=nonval, field=field, mute=True
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
            f" {colors.italic(f'{registry.__name__}()')} and save to your registry"
        )

    return result_db


def _validate(
    cls,
    values: ListLike,
    field: str | StrField | None = None,
    *,
    mute: bool = False,
    using_key: str | None = None,
    organism: str | Record | None = None,
    source: Record | None = None,
) -> np.ndarray:
    """{}"""  # noqa: D415
    from lamin_utils._inspect import validate

    return_str = True if isinstance(values, str) else False
    if isinstance(values, str):
        values = [values]
    values = _concat_lists(values)

    field = get_name_field(cls, field=field)

    queryset = _queryset(cls, using_key)
    using_key = queryset.db
    if isinstance(source, Record):
        _check_source_db(source, using_key)
        queryset = queryset.filter(source=source).all()
    _check_organism_db(organism, using_key)
    field_values = pd.Series(
        _filter_query_based_on_organism(
            queryset=queryset,
            field=field,
            organism=organism,
            values_list_field=field,
        ),
        dtype="object",
    )
    if field_values.empty:
        if not mute:
            msg = (
                f"Your {cls.__name__} registry is empty, consider populating it first!"
            )
            if hasattr(cls, "source_id"):
                msg += "\n   → use `.import_source()` to import records from a source, e.g. a public ontology"
            logger.warning(msg)
        return np.array([False] * len(values))

    result = validate(
        identifiers=values,
        field_values=field_values,
        case_sensitive=True,
        mute=mute,
        field=field,
    )
    if return_str and len(result) == 1:
        return result[0]
    else:
        return result


@classmethod  # type: ignore
@doc_args(CanValidate.standardize.__doc__)
def standardize(
    cls,
    values: ListLike,
    field: str | StrField | None = None,
    *,
    return_field: str = None,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    mute: bool = False,
    public_aware: bool = True,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    organism: str | Record | None = None,
    source: Record | None = None,
) -> list[str] | dict[str, str]:
    """{}"""  # noqa: D415
    return _standardize(
        cls=cls,
        values=values,
        field=field,
        return_field=return_field,
        return_mapper=return_mapper,
        case_sensitive=case_sensitive,
        mute=mute,
        public_aware=public_aware,
        keep=keep,
        synonyms_field=synonyms_field,
        organism=organism,
        source=source,
    )


def set_abbr(self, value: str):
    self.abbr = value

    if hasattr(self, "name") and value == self.name:
        pass
    else:
        try:
            self.add_synonym(value, save=False)
        except Exception as e:  # pragma: no cover
            logger.debug(
                f"Encountered an Exception while attempting to add synonyms.\n{e}"
            )

    if not self._state.adding:
        self.save()


def add_synonym(
    self,
    synonym: str | ListLike,
    force: bool = False,
    save: bool | None = None,
):
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(
        synonym=synonym, record=self, force=force, action="add", save=save
    )


def remove_synonym(self, synonym: str | ListLike):
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(synonym=synonym, record=self, action="remove")


def _standardize(
    cls,
    values: ListLike,
    field: str | StrField | None = None,
    *,
    return_field: str = None,
    return_mapper: bool = False,
    case_sensitive: bool = False,
    mute: bool = False,
    public_aware: bool = True,
    keep: Literal["first", "last", False] = "first",
    synonyms_field: str = "synonyms",
    using_key: str | None = None,
    organism: str | Record | None = None,
    source: Record | None = None,
) -> list[str] | dict[str, str]:
    """{}"""  # noqa: D415
    from lamin_utils._standardize import standardize as map_synonyms

    return_str = True if isinstance(values, str) else False
    if isinstance(values, str):
        values = [values]
    values = _concat_lists(values)

    field = get_name_field(cls, field=field)
    return_field = get_name_field(
        cls, field=field if return_field is None else return_field
    )
    queryset = _queryset(cls, using_key)
    using_key = queryset.db
    if isinstance(source, Record):
        _check_source_db(source, using_key)
        queryset = queryset.filter(source=source).all()
    _check_organism_db(organism, using_key)
    registry = queryset.model

    if _has_organism_field(registry):
        # here, we can safely import bionty
        from bionty._bionty import create_or_get_organism_record

        organism_record = create_or_get_organism_record(
            organism=organism, registry=registry, field=field
        )
        organism = (
            organism_record.name if organism_record is not None else organism_record
        )

    # only perform synonym mapping if field is the name field
    if hasattr(registry, "_name_field") and field != registry._name_field:
        synonyms_field = None

    try:
        registry._meta.get_field(synonyms_field)
        fields = {i for i in [field, return_field, synonyms_field] if i is not None}
        df = _filter_query_based_on_organism(
            queryset=queryset,
            field=field,
            organism=organism,
            fields=list(fields),
        )
    except FieldDoesNotExist:
        df = pd.DataFrame()

    _kwargs = {
        "field": field,
        "return_field": return_field,
        "case_sensitive": case_sensitive,
        "keep": keep,
        "synonyms_field": synonyms_field,
    }
    # standardized names from the DB
    std_names_db = map_synonyms(
        df=df,
        identifiers=values,
        return_mapper=return_mapper,
        mute=mute,
        **_kwargs,
    )

    def _return(result: list, mapper: dict):
        if return_mapper:
            return mapper
        else:
            if return_str and len(result) == 1:
                return result[0]
            return result

    # map synonyms in Bionty
    if registry.__get_schema_name__() == "bionty" and public_aware:
        mapper = {}
        if return_mapper:
            mapper = std_names_db
            std_names_db = map_synonyms(
                df=df, identifiers=values, return_mapper=False, mute=True, **_kwargs
            )

        val_res = registry.validate(
            std_names_db, field=field, mute=True, organism=organism
        )
        if all(val_res):
            return _return(result=std_names_db, mapper=mapper)

        nonval = np.array(std_names_db)[~val_res]
        std_names_bt_mapper = registry.public(organism=organism).standardize(
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
                f"\n   please add corresponding {registry._meta.model.__name__} records via"
                f" `.from_values({list(set(std_names_bt_mapper.values()))})`"
            )
            logger.warning(warn_msg)

        mapper.update(std_names_bt_mapper)
        if pd.api.types.is_categorical_dtype(std_names_db):
            result = std_names_db.cat.rename_categories(std_names_bt_mapper).tolist()
        else:
            result = pd.Series(std_names_db).replace(std_names_bt_mapper).tolist()
        return _return(result=result, mapper=mapper)

    else:
        return _return(result=std_names_db, mapper=std_names_db)


def _add_or_remove_synonyms(
    synonym: str | ListLike,
    record: Record,
    action: Literal["add", "remove"],
    force: bool = False,
    save: bool | None = None,
):
    """Add or remove synonyms."""

    def check_synonyms_in_all_records(synonyms: set[str], record: Record):
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
            raise ValidationError(
                f"you are trying to assign a synonym to record: {record}\n"
                "    → consider removing the synonym from existing records or using a different synonym."
            )

    # passed synonyms
    # nothing happens when passing an empty string or list
    if isinstance(synonym, str):
        if len(synonym) == 0:
            return
        syn_new_set = {synonym}
    else:
        if synonym == [""]:
            return
        syn_new_set = set(synonym)
    # nothing happens when passing an empty string or list
    if len(syn_new_set) == 0:
        return
    # because we use | as the separator
    if any("|" in i for i in syn_new_set):
        raise ValidationError("a synonym can't contain '|'!")

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


def _check_synonyms_field_exist(record: Record):
    try:
        record.__getattribute__("synonyms")
    except AttributeError:
        raise NotImplementedError(
            f"No synonyms field found in table {record.__class__.__name__}!"
        ) from None


def _filter_query_based_on_organism(
    queryset: QuerySet,
    field: str,
    organism: str | Record | None = None,
    values_list_field: str | None = None,
    fields: list[str] | None = None,
):
    """Filter a queryset based on organism."""
    import pandas as pd

    registry = queryset.model

    if _has_organism_field(registry) and not _field_is_id(field, registry):
        # here, we can safely import bionty
        from bionty._bionty import create_or_get_organism_record

        organism_record = create_or_get_organism_record(
            organism=organism, registry=registry, field=field
        )
        if organism_record is not None:
            queryset = queryset.filter(organism__name=organism_record.name)

    if values_list_field is None:
        if fields:
            return pd.DataFrame.from_records(
                queryset.values_list(*fields), columns=fields
            )
        return pd.DataFrame.from_records(queryset.values())

    else:
        return queryset.values_list(values_list_field, flat=True)


def _field_is_id(field: str, registry: type[Record]) -> bool:
    """Check if the field is an ontology ID."""
    if hasattr(registry, "_ontology_id_field"):
        if field == registry._ontology_id_field:
            return True
    if field.endswith("id"):
        return True
    return False


METHOD_NAMES = [
    "validate",
    "inspect",
    "standardize",
    "add_synonym",
    "remove_synonym",
    "set_abbr",
    "from_values",
]

if ln_setup._TESTING:  # type: ignore
    from inspect import signature

    SIGS = {
        name: signature(getattr(CanValidate, name))
        for name in METHOD_NAMES
        if not name.startswith("__")
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, CanValidate, globals())
