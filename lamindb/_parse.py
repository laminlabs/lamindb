from typing import Any, Dict, Iterable, List, Optional, Tuple, TypeVar, Union

import numpy as np
import pandas as pd
from django.db.models import Model
from django.db.models.query_utils import DeferredAttribute as Field
from lamin_logger import logger

from ._select import select

ListLike = TypeVar("ListLike", pd.Series, list, np.array)


def parse(
    iterable: Union[ListLike, pd.DataFrame],
    field: Union[Field, Dict[str, Field]],
    species: str = None,
) -> List[Model]:
    """Parse a dataset column based on a Model entity field.

    Guide: :doc:`/biology/ontologies`.

    Args:
        iterable: a `ListLike` of values or a `DataFrame`.
        field: if iterable is `ListLike`: a `Model` field to parse into.
               if iterable is `DataFrame`: a dict of {column_name : Model_field}.
        species: if None, will use default species in bionty for each entity.

    Returns:
        A list of Model records.
    """
    if isinstance(iterable, pd.DataFrame):
        logger.warning(
            "Providing a DataFrame will always lead to creating new records!"
        )
        if not isinstance(field, dict):
            raise TypeError("field must be a dictionary of {column_name: Field}!")
        column_mapper = {colname: f.field.name for colname, f in field.items()}
        class_mapper = {f.field.name: f.field.model for f in field.values()}
        if len(set(class_mapper.values())) > 1:
            raise NotImplementedError("fields must from the same entity!")
        model = list(class_mapper.values())[0]
        df = iterable.rename(columns=column_mapper)
        df = df.dropna().drop_duplicates()
        df = df.mask(df == "", None)
        model_field_names = {i.name for i in model._meta.fields}
        df = df.loc[:, df.columns.isin(model_field_names)]
        if df.index.name is not None:
            df = df.reset_index()
        iterable = df.to_dict(orient="records")
        return [model(**kwargs) for kwargs in iterable]
    else:
        if not isinstance(field, Field):
            raise TypeError("field must be an ORM field, e.g., `CellType.name`!")
        iterable = set(iterable)
        return get_or_create_records(iterable=iterable, field=field, species=species)


def not_null_iterable(iterable: Iterable) -> set:
    # No entries are made for NAs, '', None
    return set([i for i in iterable if (not pd.isnull(i)) and (i != "")])


def get_or_create_records(
    iterable: Iterable,
    field: Field,
    species: Optional[str] = None,
) -> List:
    """Get or create records from iterables."""
    iterable = not_null_iterable(iterable)
    model = field.field.model  # is DeferredAttribute
    field_name = field.field.name

    kwargs_species, condition_species, bionty_df = _preprocess_species(
        species=species, model=model
    )

    # get all existing records in the db for a species
    condition = {f"{field_name}__in": iterable}
    condition.update(condition_species)
    stmt = select(model, **condition)

    records = stmt.list()
    if len(records) > 0:
        logger.hint(f"Returned {len(records)} existing {model.__name__} records.")

    existing_values = iterable.intersection(stmt.values_list(field_name, flat=True))

    # new records to be created based on new values
    new_values = iterable.difference(existing_values)
    if len(new_values) > 0:
        # first try to populate additional fields from bionty
        mapped_values = bionty_df[field_name].intersection(new_values)
        if len(mapped_values) > 0:
            new_values_kwargs = _bulk_create_dicts_from_df(
                keys=mapped_values, column_name=field_name, df=bionty_df
            )
            for kwargs in new_values_kwargs:
                kwargs.update(kwargs_species)
                records.append(model(**kwargs))
            logger.hint(
                f"Created {len(mapped_values)} {model.__name__} records from bionty"
                f" with {bionty_df.shape[1]} fields"
            )
        # unmapped new_ids will only create records with field and species
        unmapped_values = new_values.difference(mapped_values)
        if len(unmapped_values) > 0:
            for i in unmapped_values:
                kwargs = {field_name: i}
                kwargs.update(kwargs_species)
                records.append(model(**kwargs))
            logger.hint(
                f"Created {len(unmapped_values)} {model.__name__} records with a single"
                f" field '{field_name}'"
            )
    return records


def _preprocess_species(
    model: Model, species: Optional[str] = None
) -> Tuple[dict, dict, pd.DataFrame]:
    kwargs_species: Dict = {}
    condition_species: Dict = {}
    from lnschema_bionty._bionty import create_or_get_species_record, get_bionty_object

    bionty_object = get_bionty_object(model=model, species=species)

    if bionty_object is not None:
        species_record = create_or_get_species_record(species=bionty_object.species)

        # if species is specified, only pull species-specific records
        condition_species["species__name"] = species_record.name
        # add species to the record
        kwargs_species["species"] = species_record
        logger.info(
            f"Returning records with species='{species_record.name}'."  # type:ignore
        )
    bionty_df = _filter_bionty_df_columns(model=model, bionty_object=bionty_object)

    return kwargs_species, condition_species, bionty_df


def _filter_bionty_df_columns(model: Model, bionty_object: Any) -> pd.DataFrame:
    bionty_df = pd.DataFrame()
    if bionty_object is not None:
        model_field_names = {i.name for i in model._meta.fields}
        bionty_df = bionty_object.df().reset_index()
        bionty_df = bionty_df.loc[:, bionty_df.columns.isin(model_field_names)]
    return bionty_df


def _bulk_create_dicts_from_df(keys: list, column_name: str, df: pd.DataFrame) -> dict:
    """Get fields from a dataframe for many rows."""
    if df.index.name != column_name:
        df = df.set_index(column_name)
    # keep the last record (assuming most recent) if duplicated
    df = df[~df.index.duplicated(keep="last")]
    return df.loc[keys].reset_index().to_dict(orient="records")
