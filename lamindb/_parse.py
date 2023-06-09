from typing import Any, Dict, Iterable, List, Optional, Tuple, TypeVar, Union

import numpy as np
import pandas as pd
from django.core.exceptions import MultipleObjectsReturned, ObjectDoesNotExist
from django.db.models import Model, Q
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_logger import colors, logger

from ._select import select

ListLike = TypeVar("ListLike", pd.Series, list, np.array)


def parse(
    iterable: Union[ListLike, pd.DataFrame],
    field: Union[Field, Dict[str, Field]],
    species: str = None,
    query_existing: bool = True,
) -> List[Model]:
    """Parse a dataset column based on a Model entity field.

    Guide: :doc:`/biology/ontologies`.

    Args:
        iterable: a `ListLike` of values or a `DataFrame`.
        field: if iterable is `ListLike`: a `Model` field to parse into.
               if iterable is `DataFrame`: a dict of {column_name : Model_field}.
        species: if None, will use default species in bionty for each entity.
        query_existing: if False, always create new records

    Returns:
        A list of Model records.
    """
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

        if not query_existing:
            # always create new records, skips query for existing
            records = [model(**kwargs) for kwargs in df_records]
            n_new = len(records)
            logger.hint(f"Created {n_new} {model.__name__} records")
            return records

        # make sure to only return 1 existing entry for each row
        queryset = _bulk_query_fields(df_records=df_records, model=model)
        n_exist, n_new, records = _get_from_queryset(
            queryset=queryset, df_records=df_records, model=model
        )

        if n_exist > 0:
            text = colors.green(f"{n_exist} existing {model.__name__} records")
            logger.hint(f"Returned {text}")
        if n_new > 0:
            text = colors.purple(f"{n_new} {model.__name__} records")
            logger.hint(f"Created {text} with {df.shape[1]} fields")
        return records
    else:
        if not isinstance(field, Field):
            raise TypeError("field must be an ORM field, e.g., `CellType.name`!")
        return get_or_create_records(iterable=iterable, field=field, species=species)


def index_iterable(iterable: Iterable) -> pd.Index:
    idx = pd.Index(iterable).unique()
    # No entries are made for NAs, '', None
    # returns an ordered unique not null list
    return idx[(idx != "") & (~idx.isnull())]


def get_or_create_records(
    iterable: Iterable,
    field: Field,
    species: Optional[str] = None,
) -> List:
    """Get or create records from iterables."""
    iterable_idx = index_iterable(iterable)
    model = field.field.model  # is DeferredAttribute
    field_name = field.field.name

    kwargs_species_source, condition_species, bionty_df = _preprocess_species(
        species=species, model=model
    )

    # get all existing records in the db for a species
    condition = {f"{field_name}__in": iterable}
    condition.update(condition_species)
    stmt = select(model, **condition)

    records = stmt.list()
    if len(records) > 0:
        text = colors.green(f"{len(records)} existing {model.__name__} records")
        logger.hint(f"Returned {text}")

    existing_values = iterable_idx.intersection(stmt.values_list(field_name, flat=True))

    # new records to be created based on new values
    new_values = iterable_idx.difference(existing_values)
    if len(new_values) > 0:
        # first try to populate additional fields from bionty
        mapped_values = new_values.intersection(bionty_df[field_name])
        if len(mapped_values) > 0:
            new_values_kwargs = _bulk_create_dicts_from_df(
                keys=mapped_values, column_name=field_name, df=bionty_df
            )
            for kwargs in new_values_kwargs:
                kwargs.update(kwargs_species_source)
                records.append(model(**kwargs))
            text = colors.purple(f"{len(mapped_values)} {model.__name__} records")
            logger.hint(f"Created {text} from bionty with {bionty_df.shape[1]} fields")
            bionty_version = kwargs_species_source.get("bionty_source")
            if bionty_version is not None:
                logger.hint(f"Using bionty source: {bionty_version}")
        # unmapped new_ids will only create records with field and species
        unmapped_values = new_values.difference(mapped_values)
        if len(unmapped_values) > 0:
            for i in unmapped_values:
                kwargs = {field_name: i}
                species_record = kwargs_species_source.get("species")
                if species_record is not None:
                    kwargs.update({"species": species_record})
                records.append(model(**kwargs))
            text = colors.blue(f"{len(unmapped_values)} {model.__name__} records")
            logger.hint(f"Created {text} records with a single field '{field_name}'")
    return records


def _get_from_queryset(queryset, df_records, model):
    n_exist = 0
    n_new = 0
    records = []
    for kwargs in df_records:
        try:
            records.append(queryset.get(**kwargs))
            n_exist += 1
        except MultipleObjectsReturned:
            records.append(queryset.filter(**kwargs).first())
            logger.warning(
                f"Found multiple existing {model.__name__} records with"
                f" {kwargs}, returning the first query result!"
            )
            n_exist += 1
        except ObjectDoesNotExist:
            records.append(model(**kwargs))
            n_new += 1
    return n_exist, n_new, records


def _preprocess_species(
    model: Model, species: Optional[str] = None
) -> Tuple[dict, dict, pd.DataFrame]:
    kwargs_species_source: Dict = {}
    condition_species: Dict = {}
    from lnschema_bionty._bionty import (
        create_or_get_species_record,
        get_bionty_object,
        get_bionty_source_record,
    )

    bionty_object = get_bionty_object(model=model, species=species)

    if bionty_object is not None:
        species_record = create_or_get_species_record(species=bionty_object.species)
        if species_record is not None:
            # if species is specified, only pull species-specific records
            condition_species["species__name"] = species_record.name
            # add species to the record
            kwargs_species_source["species"] = species_record
            logger.info(f"Returning records with species='{species_record.name}'...")
        kwargs_species_source["bionty_source"] = get_bionty_source_record(bionty_object)
    bionty_df = _filter_bionty_df_columns(model=model, bionty_object=bionty_object)

    return kwargs_species_source, condition_species, bionty_df


def _filter_bionty_df_columns(model: Model, bionty_object: Any) -> pd.DataFrame:
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


def _bulk_query_fields(df_records: list, model: Model):
    condition = Q(**df_records[0])
    for kwargs in df_records[1:]:
        condition = condition.__getattribute__("__or__")(Q(**kwargs))

    queryset = model.objects.filter(condition)

    return queryset


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
