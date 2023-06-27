from typing import Any, Dict, Iterable, List, Optional, TypeVar, Union

import numpy as np
import pandas as pd
from django.db.models import Q
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_logger import colors, logger
from lamindb_setup.dev import deprecated
from lnschema_core.models import BaseORM

from ._select import select
from .dev._settings import settings

ListLike = TypeVar("ListLike", pd.Series, list, np.array)


# The base function for `from_iter` and `from_bionty`
def get_or_create_records(
    iterable: ListLike,
    field: Field,
    *,
    from_bionty: bool = False,
    **kwargs,
) -> List:
    """Get or create records from iterables."""
    upon_create_search_names = settings.upon_create_search_names
    settings.upon_create_search_names = False
    try:
        field_name = field.field.name
        model = field.field.model
        iterable_idx = index_iterable(iterable)

        records, nonexist_values = get_existing_records(
            iterable_idx=iterable_idx, field=field, kwargs=kwargs
        )

        # new records to be created based on new values
        if len(nonexist_values) > 0:
            if from_bionty:
                records_bionty, unmapped_values = create_records_from_bionty(
                    iterable_idx=nonexist_values, field=field, **kwargs
                )
                records += records_bionty
            else:
                unmapped_values = nonexist_values
            # unmapped new_ids will only create records with field and kwargs
            if len(unmapped_values) > 0:
                for i in unmapped_values:
                    records.append(model(**{field_name: i}, **kwargs))
                logger.hint(
                    "Created"
                    f" {colors.red(f'{len(unmapped_values)} {model.__name__} records')}"
                    f" with a single field {colors.red(f'{field_name}')}"
                )
        return records
    finally:
        settings.upon_create_search_names = upon_create_search_names


@deprecated("ORM.from_iter()")
def parse(
    iterable: Union[ListLike, pd.DataFrame],
    field: Union[Field, Dict[str, Field]],
    *,
    species: Optional[str] = None,
) -> List[BaseORM]:
    """Parse identifiers and create records through lookups for a given field.

    Guide: :doc:`/biology/registries`.

    Args:
        iterable: `Union[ListLike, pd.DataFrame]` A `ListLike` of identifiers or
            a `DataFrame`.
        field: `Union[Field, Dict[str, Field]]` If `iterable` is `ListLike`, a
            `BaseORM` field to look up.
            If `iterable` is `DataFrame`, a dict of `{column_name1: field1,
            column_name2: field2}`.
        species: `Optional[str]` Either `"human"`, `"mouse"`, or any other
            `name` of `Bionty.Species`. If `None`, will use default species in
            bionty for each entity.

    Returns:
        A list of records.

    For every `value` in an iterable of identifiers and a given `ORM.field`,
    this function performs:

    1. It checks whether the value already exists in the database
       (`ORM.select(field=value)`). If so, it adds the queried record to
       the returned list and skips step 2. Otherwise, proceed with 2.
    2. If the `ORM` is from `lnschema_bionty`, it checks whether there is an
       exact match in the underlying ontology (`Bionty.inspect(value, field)`).
       If so, it creates a record from Bionty and adds it to the returned list.
       Otherwise, it create a record that populates a single field using `value`
       and adds the record to the returned list.

    """
    upon_create_search_names = settings.upon_create_search_names
    settings.upon_create_search_names = False
    try:
        if isinstance(iterable, pd.DataFrame):
            # check the field must be a dictionary
            if not isinstance(field, dict):
                raise TypeError("field must be a dictionary of {column_name: Field}!")

            # check only one single model class is passed
            class_mapper = {f.field.name: f.field.model for f in field.values()}
            if len(set(class_mapper.values())) > 1:
                raise NotImplementedError("fields must from the same entity!")
            model = list(class_mapper.values())[0]

            df = _map_columns_to_fields(df=iterable, field=field)
            df_records = df.to_dict(orient="records")

            # make sure to only return 1 existing entry for each row
            queryset = get_existing_records_multifields(
                df_records=df_records, model=model
            )
            records = queryset.list()
            df_records_new = [
                i for i in df_records if not queryset.filter(**i).exists()
            ]

            if len(records) > 0:
                logger.hint(
                    "Returned"
                    f" {colors.green(f'{len(records)} existing {model.__name__} DB records')}"  # noqa
                )
            if len(df_records_new) > 0:
                logger.hint(
                    "Created"
                    f" {colors.purple(f'{len(df_records_new)} {model.__name__} records')} with"  # noqa
                    f" {df.shape[1]} fields"
                )
                records += [model(**i) for i in df_records_new]
            return records
        else:
            if not isinstance(field, Field):
                raise TypeError("field must be an ORM field, e.g., `CellType.name`!")
            return get_or_create_records(
                iterable=iterable, field=field, species=species
            )
    finally:
        settings.upon_create_search_names = upon_create_search_names


def index_iterable(iterable: Iterable) -> pd.Index:
    idx = pd.Index(iterable).unique()
    # No entries are made for NAs, '', None
    # returns an ordered unique not null list
    return idx[(idx != "") & (~idx.isnull())]


def get_existing_records(iterable_idx: pd.Index, field: Field, kwargs: Dict = {}):
    field_name = field.field.name
    model = field.field.model

    # map synonyms based on the DB reference
    syn_mapper = model.map_synonyms(
        iterable_idx, species=kwargs.get("species"), return_mapper=True
    )

    syn_msg = ""
    if len(syn_mapper) > 0:
        syn_msg = (
            "Returned"
            f" {colors.green(f'{len(syn_mapper)} existing {model.__name__} DB records')} that"  # noqa
            f" matched {colors.green('synonyms')}"
        )
        iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index

    # get all existing records in the db
    # if necessary, create records for the values in kwargs
    # k:v -> k:v_record
    # kwargs is used to deal with species
    condition = {f"{field_name}__in": iterable_idx.values}
    kwargs, condition = _species_kwargs(orm=model, kwargs=kwargs, condition=condition)

    stmt = select(model, **condition)

    records = stmt.list()  # existing records
    n_name = len(records) - len(syn_mapper)
    if n_name > 0:
        logger.hint(
            "Returned"
            f" {colors.green(f'{n_name} existing {model.__name__} DB records')} that"
            f" matched {colors.green(f'{field_name}')} field"
        )
    # make sure that synonyms logging appears after the field logging
    if len(syn_msg) > 0:
        logger.hint(syn_msg)

    existing_values = iterable_idx.intersection(stmt.values_list(field_name, flat=True))
    nonexist_values = iterable_idx.difference(existing_values)

    return records, nonexist_values


def get_existing_records_multifields(
    df_records: List, model: BaseORM, kwargs: Dict = {}
):
    q = Q(**df_records[0])
    for df_record in df_records[1:]:
        q = q.__getattribute__("__or__")(Q(**df_record))

    kwargs, condition = _species_kwargs(orm=model, kwargs=kwargs)
    stmt = model.select(**condition).filter(q)
    return stmt


def _species_kwargs(orm: BaseORM, kwargs: Dict = {}, condition: Dict = {}):
    """Create records based on the kwargs."""
    if kwargs.get("species") is not None:
        from lnschema_bionty._bionty import create_or_get_species_record

        species_record = create_or_get_species_record(
            species=kwargs.get("species"), orm=orm
        )
        if species_record is not None:
            kwargs.update({"species": species_record})
            condition.update({"species__name": species_record.name})

    return kwargs, condition


def create_records_from_bionty(
    iterable_idx: pd.Index,
    field: Field,
    **kwargs,
):
    model = field.field.model
    field_name = field.field.name
    records: List = []
    # populate additional fields from bionty
    from lnschema_bionty._bionty import get_bionty_object, get_bionty_source_record

    # create the corresponding bionty object from model
    bionty_object = get_bionty_object(orm=model, species=kwargs.get("species"))
    # add bionty_source record to the kwargs
    kwargs.update({"bionty_source": get_bionty_source_record(bionty_object)})

    # filter the columns in bionty df based on fields
    bionty_df = _filter_bionty_df_columns(model=model, bionty_object=bionty_object)

    # map synonyms in the bionty reference
    try:
        syn_mapper = bionty_object.map_synonyms(iterable_idx, return_mapper=True)
    except KeyError:
        # no synonyms column
        syn_mapper = {}
    msg_syn: str = ""
    if len(syn_mapper) > 0:
        msg_syn = (
            "Created"
            f" {colors.purple(f'{len(syn_mapper)} {model.__name__} records from Bionty')} that"  # noqa
            f" matched {colors.purple('synonyms')}"
        )
        iterable_idx = iterable_idx.to_frame().rename(index=syn_mapper).index

    # create records for values that are found in the bionty reference
    mapped_values = iterable_idx.intersection(bionty_df[field_name])

    if len(mapped_values) > 0:
        bionty_kwargs = _bulk_create_dicts_from_df(
            keys=mapped_values, column_name=field_name, df=bionty_df
        )
        for bk in bionty_kwargs:
            records.append(model(**bk, **kwargs))

        # logging of BiontySource linking
        source_msg = (
            ""
            if kwargs.get("bionty_source") is None
            else f", linked to BiontySource id={kwargs.get('bionty_source').id}"  # type:ignore # noqa
        )

        # number of records that matches field (not synonyms)
        n_name = len(records) - len(syn_mapper)
        if n_name > 0:
            msg = (
                "Created"
                f" {colors.purple(f'{n_name} {model.__name__} records from Bionty')} that"  # noqa
                f" matched {colors.purple(f'{field_name}')} field"
            )
            logger.hint(msg + source_msg)
        # make sure that synonyms logging appears after the field logging
        if len(msg_syn) > 0:
            logger.hint(msg_syn + source_msg)

    # return the values that are not found in the bionty reference
    unmapped_values = iterable_idx.difference(mapped_values)
    return records, unmapped_values


def _filter_bionty_df_columns(model: BaseORM, bionty_object: Any) -> pd.DataFrame:
    bionty_df = pd.DataFrame()
    if bionty_object is not None:
        model_field_names = {i.name for i in model._meta.fields}
        bionty_df = bionty_object.df().reset_index()
        bionty_df = bionty_df.loc[:, bionty_df.columns.isin(model_field_names)]
    return bionty_df


def _bulk_create_dicts_from_df(
    keys: Union[set, List], column_name: str, df: pd.DataFrame
) -> dict:
    """Get fields from a DataFrame for many rows."""
    if df.index.name != column_name:
        df = df.set_index(column_name)
    # keep the last record (assuming most recent) if duplicated
    df = df[~df.index.duplicated(keep="last")]
    return df.loc[list(keys)].reset_index().to_dict(orient="records")


def _map_columns_to_fields(df: pd.DataFrame, field: dict) -> pd.DataFrame:
    """Subset dataframe to mappable fields columns and clean up."""
    column_mapper = {colname: f.field.name for colname, f in field.items()}
    # subset to columns containing fields
    df = df.copy()
    if df.index.name is not None:
        df = df.reset_index()
    df = df.loc[:, df.columns.isin(field.keys())]
    df = df.rename(columns=column_mapper)
    df = df.dropna().drop_duplicates()
    # TODO: remove after having the auto conversion for django ORMs
    df = df.mask(df == "", None)
    return df
