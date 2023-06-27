import builtins
from typing import (
    Dict,
    Iterable,
    List,
    Literal,
    NamedTuple,
    Optional,
    Set,
    Tuple,
    Union,
)

import pandas as pd
from django.core.exceptions import FieldDoesNotExist
from django.db.models import CharField, TextField
from lamin_logger import logger
from lamin_logger._lookup import Lookup
from lnschema_core import BaseORM

from ._from_values import Field, ListLike, get_or_create_records
from .dev._settings import settings

_is_ipython = getattr(builtins, "__IPYTHON__", False)


class ValidationError(Exception):
    pass


def validate_required_fields(orm: BaseORM, kwargs):
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


def suggest_objects_with_same_name(orm: BaseORM, kwargs) -> Optional[str]:
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


def return_object_from_bionty(orm: BaseORM, *args, **kwargs) -> Dict:
    """Pass bionty search/lookup results."""
    from lnschema_bionty._bionty import (
        create_or_get_species_record,
        encode_id,
        get_bionty_source_record,
    )

    arg = args[0]
    if isinstance(arg, Tuple):  # type:ignore
        bionty_kwargs = arg._asdict()
    else:
        bionty_kwargs = arg[0]._asdict()

    if len(bionty_kwargs) > 0:
        import bionty as bt

        # add species and bionty_source
        species_record = create_or_get_species_record(
            orm=orm, species=kwargs.get("species")
        )
        if species_record is not None:
            bionty_kwargs["species"] = species_record
        bionty_object = getattr(bt, orm.__class__.__name__)(
            species=species_record.name if species_record is not None else None
        )
        bionty_kwargs["bionty_source"] = get_bionty_source_record(bionty_object)

        model_field_names = {i.name for i in orm._meta.fields}
        bionty_kwargs = {
            k: v for k, v in bionty_kwargs.items() if k in model_field_names
        }
    return encode_id(orm=orm, kwargs=bionty_kwargs)


def __init__(orm: BaseORM, *args, **kwargs):
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
                super(BaseORM, orm).__init__(*new_args)
                orm._state.adding = False  # mimic from_db
                return None
        if orm.__module__.startswith("lnschema_bionty"):
            from lnschema_bionty._bionty import encode_id

            kwargs = encode_id(orm=orm, kwargs=kwargs)
        super(BaseORM, orm).__init__(**kwargs)
    elif (
        orm.__module__.startswith("lnschema_bionty")
        and args
        and len(args) == 1
        and isinstance(args[0], (Tuple, List))  # type:ignore
        and len(args[0]) > 0
    ):
        if isinstance(args[0], List) and len(args[0]) > 1:
            logger.warning(
                "Multiple lookup/search results are passed, only returning record from"
                " the first entry"
            )
        result = return_object_from_bionty(orm, *args, **kwargs)  # type:ignore
        try:
            existing_object = orm.select(**result)[0]
            new_args = [
                getattr(existing_object, field.attname)
                for field in orm._meta.concrete_fields
            ]
            super(BaseORM, orm).__init__(*new_args)
            orm._state.adding = False  # mimic from_db
        except IndexError:
            super(BaseORM, orm).__init__(**result)
    elif len(args) != len(orm._meta.concrete_fields):
        raise ValueError("Please provide keyword arguments, not plain arguments")
    else:
        # object is loaded from DB (**kwargs could be ommitted below, I believe)
        super(BaseORM, orm).__init__(*args, **kwargs)


@classmethod  # type:ignore
def from_values(cls, values: ListLike, field: Union[Field, str], **kwargs):
    if isinstance(field, str):
        field = getattr(cls, field)
    if not isinstance(field, Field):  # field is DeferredAttribute
        raise TypeError(
            "field must be a string or an ORM field, e.g., `CellType.name`!"
        )
    if cls.__name__ == "FeatureSet":
        from lamindb._featureset_methods import parse_features_from_iterable

        features = parse_features_from_iterable(
            iterable=values,
            field=field,
            species=kwargs.get("species"),
        )
        return features

    from_bionty = True if cls.__module__.startswith("lnschema_bionty.") else False
    return get_or_create_records(
        iterable=values, field=field, from_bionty=from_bionty, **kwargs
    )


@classmethod  # type: ignore
def search(
    cls,
    string: str,
    *,
    field: Optional[Union[str, CharField, TextField]] = None,
    top_hit: bool = False,
    case_sensitive: bool = True,
    synonyms_field: Optional[Union[str, TextField, CharField]] = "synonyms",
    synonyms_sep: str = "|",
) -> Union[pd.DataFrame, BaseORM]:
    """Search the table.

    Args:
        string: `str` The input string to match against the field ontology values.
        field: `Optional[Union[str, CharField, TextField]] = None` The field
            against which the input string is matching.
        top_hit: `bool = False` If `True`, return only the top hit or hits (in
            case of equal scores).
        case_sensitive: `bool = False` Whether the match is case sensitive.
        synonyms_field: `bool = True` Also search synonyms. If `None`, is ignored.

    Returns:
        A sorted `DataFrame` of search results with a score in column
        `__ratio__`. If `top_hit` is `True`, the best match.
    """
    import pandas as pd
    from lamin_logger._search import search

    if field is None:
        field = get_default_str_field(cls)
    if not isinstance(field, str):
        field = field.field.name

    records = cls.objects.all()
    df = pd.DataFrame.from_records(records.values())

    result = search(
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
def lookup(cls, field: Optional[Union[str, CharField, TextField]] = None) -> NamedTuple:
    """Return an auto-complete object for a field.

    Args:
        field: `Optional[Union[str, CharField, TextField]] = None` The field to
            look up the values for. Defaults to 'name'.

    Returns:
        A `NamedTuple` of lookup information of the field values with a
        dictionary converter.

    Examples:
        >>> import lnschema_bionty as lb
        >>> lookup = lb.Gene.lookup()
        >>> lookup.adgb_dt
        >>> lookup_dict = lookup.dict()
        >>> lookup['ADGB-DT']
    """
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


lookup.__doc__ = Lookup.__doc__


@classmethod  # type: ignore
def inspect(
    cls,
    identifiers: Iterable,
    field: Union[str, CharField, TextField],
    *,
    case_sensitive: bool = False,
    inspect_synonyms: bool = True,
    return_df: bool = False,
    logging: bool = True,
    **kwargs,
) -> Union[pd.DataFrame, Dict[str, List[str]]]:
    """Inspect if a list of identifiers are mappable to existing values of a field.

    Args:
        identifiers: Identifiers that will be checked against the field.
        field: `Union[str, CharField, TextField]` The field of identifiers.
                Examples are 'ontology_id' to map against the source ID
                or 'name' to map against the ontologies field names.
        case_sensitive: Whether the identifier inspection is case sensitive.
        inspect_synonyms: Whether to inspect synonyms.
        return_df: Whether to return a Pandas DataFrame.

    Returns:
        - A Dictionary of "mapped" and "unmapped" identifiers
        - If `return_df`: A DataFrame indexed by identifiers with a boolean `__mapped__`
            column that indicates compliance with the identifiers.

    Examples:
        >>> import lnschema_bionty as lb
        >>> gene_symbols = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
        >>> lb.Gene.inspect(gene_symbols, field=lb.Gene.symbol)
    """
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
    """Maps input synonyms to standardized names.

    Args:
        synonyms: `Iterable` Synonyms that will be standardized.
        return_mapper: `bool = False` If `True`, returns `{input_synonym1:
            standardized_name1}`.
        case_sensitive: `bool = False` Whether the mapping is case sensitive.
        species: `Optional[str]` Map only against this species related entries.
        keep: `Literal["first", "last", False] = "first"` When a synonym maps to
            multiple names, determines which duplicates to mark as
            `pd.DataFrame.duplicated`

                - "first": returns the first mapped standardized name
                - "last": returns the last mapped standardized name
                - `False`: returns all mapped standardized name
        synonyms_field: `str = "synonyms"` A field containing the concatenated synonyms.
        synonyms_sep: `str = "|"` Which separator is used to separate synonyms.
        field: `Optional[str]` The field representing the standardized names.

    Returns:
        If `return_mapper` is `False`: a list of standardized names. Otherwise,
        a dictionary of mapped values with mappable synonyms as keys and
        standardized names as values.

    Examples:
        >>> import lnschema_bionty as lb
        >>> gene_synonyms = ["A1CF", "A1BG", "FANCD1", "FANCD20"]
        >>> standardized_names = lb.Gene.map_synonyms(gene_synonyms, species="human")
    """
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


def _filter_df_based_on_species(
    orm: BaseORM, species: Union[str, BaseORM, None] = None
):
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
        elif isinstance(species, BaseORM):
            species_name = species.name
        else:
            species_name = species
        records = records.filter(species__name=species_name)
    except FieldDoesNotExist:
        pass

    return pd.DataFrame.from_records(records.values())


def get_default_str_field(orm: BaseORM) -> str:
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
    record: BaseORM,
    action: Literal["add", "remove"],
    force: bool = False,
):
    """Add or remove synonyms."""

    def check_synonyms_in_all_records(synonyms: Set[str], record: BaseORM):
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


def _check_synonyms_field_exist(record: BaseORM):
    try:
        record.__getattribute__("synonyms")
    except AttributeError:
        raise NotImplementedError(
            f"No synonyms field found in table {record.__class__.__name__}!"
        )


def add_synonym(self, synonym: Union[str, Iterable], force: bool = False):
    """Add synonyms to a record."""
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(synonym=synonym, record=self, force=force, action="add")


def remove_synonym(self, synonym: Union[str, Iterable]):
    """Remove synonyms from a record."""
    _check_synonyms_field_exist(self)
    _add_or_remove_synonyms(synonym=synonym, record=self, action="remove")


BaseORM.__init__ = __init__
BaseORM.search = search
BaseORM.lookup = lookup
BaseORM.map_synonyms = map_synonyms
BaseORM.inspect = inspect
BaseORM.add_synonym = add_synonym
BaseORM.remove_synonym = remove_synonym
BaseORM.from_values = from_values
