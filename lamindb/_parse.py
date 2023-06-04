from typing import Dict, Iterable, List, TypeVar, Union

import numpy as np
import pandas as pd
from lamin_logger import logger
from sqlalchemy.orm.attributes import InstrumentedAttribute
from sqlmodel import SQLModel

from ._add import add
from ._select import select

ListLike = TypeVar("ListLike", pd.Series, list, np.array)


def parse(
    iterable: Union[ListLike, pd.DataFrame],
    field: Union[InstrumentedAttribute, Dict[str, InstrumentedAttribute]],
    species: str = None,
) -> List[SQLModel]:
    """Parse a dataset column based on a SQLModel entity field.

    Guide: :doc:`/biology/ontologies`.

    Args:
        iterable: a `ListLike` of values or a `DataFrame`.
        field: if iterable is `ListLike`: a `SQLModel` field to parse into.
               if iterable is `DataFrame`: a dict of {column_name : SQLModel_field}.
        species: if None, will use default species in bionty for each entity.

    Returns:
        A list of SQLModel records.
    """
    if isinstance(iterable, pd.DataFrame):
        logger.warning(
            "Providing a DataFrame will always lead to creating new records!"
        )
        if not isinstance(field, dict):
            raise TypeError(
                "field must be a dictionary of {column_name : SQLModel_field}!"
            )
        column_mapper = {colname: f.name for colname, f in field.items()}
        class_mapper = {f.name: f.class_ for f in field.values()}
        if len(set(class_mapper.values())) > 1:
            raise NotImplementedError("fields must from the same entity!")
        model = list(class_mapper.values())[0]
        df = iterable.rename(columns=column_mapper)
        df = df.dropna().drop_duplicates()
        df = df.mask(df == "", None)
        iterable = df.reset_index().to_dict(orient="records")
        return [model(**kwargs) for kwargs in iterable]
    else:
        if not isinstance(field, InstrumentedAttribute):
            raise TypeError("field must be a SQLModel field!")
        iterable = set(iterable)
        return get_or_create_records(iterable=iterable, field=field, species=species)


def get_or_create_records(
    iterable: Iterable,
    field: InstrumentedAttribute,
    species: str = None,
) -> List:
    # No entries are made for NAs, '', None
    iterable = set([i for i in iterable if (not pd.isnull(i)) and (i != "")])
    model = field.class_
    parsing_id = field.name

    # get all existing records in the db
    # if species is specified, only pull species-specific records
    stmt = select(model).where(getattr(model, parsing_id).in_(iterable))

    # for bionty records, will add species if needed
    additional_kwargs = {}
    reference_df = pd.DataFrame()
    if model.__module__.startswith("lnschema_bionty."):
        import bionty as bt

        if species is None or species == "all":
            bionty_object = getattr(bt, model.__name__)()
        else:
            bionty_object = getattr(bt, model.__name__)(species=species)

        # insert species entry if not exists
        if bionty_object.species != "all":
            from lnschema_bionty import Species

            species = select(Species, name=bionty_object.species).one_or_none()
            if species is None:
                species = add(Species.from_bionty(name=bionty_object.species))
            try:
                stmt = stmt.where(
                    getattr(model, "species_id") == species.id  # type:ignore
                )
                additional_kwargs = {"species_id": species.id}  # type:ignore
                logger.info(
                    f"Returned records with species='{species.name}'."  # type:ignore
                )
            except AttributeError:
                pass

        reference_df = bionty_object.df().reset_index().set_index(parsing_id)

    records = stmt.all()
    if len(records) > 0:
        logger.hint(f"Returned {len(records)} existing records.")

    existing_ids = iterable.intersection(stmt.df()[parsing_id])

    # new records to be appended
    new_ids = iterable.difference(existing_ids)
    if len(new_ids) > 0:
        # mapped new_ids will fetch fields values from bionty
        mapped_id = reference_df.index.intersection(new_ids)
        if len(mapped_id) > 0:
            logger.hint(f"Created {len(mapped_id)} records from bionty.")
            new_ids_kwargs = _bulk_create_dicts_from_df(
                keys=mapped_id, column_name=parsing_id, df=reference_df
            )
            for kwargs in new_ids_kwargs:
                kwargs.update(additional_kwargs)
                records.append(model(**kwargs))
        # unmapped new_ids will only create records with parsing_id and species
        unmapped = set(new_ids).difference(mapped_id)
        if len(unmapped) > 0:
            logger.hint(f"Created {len(unmapped)} records from user inputs.")
            for i in unmapped:
                kwargs = {parsing_id: i}
                kwargs.update(additional_kwargs)
                records.append(model(**kwargs))
    return records


def _bulk_create_dicts_from_df(keys: list, column_name: str, df: pd.DataFrame) -> dict:
    """Get fields from a dataframe for many rows."""
    if df.index.name != column_name:
        df = df.set_index(column_name)
    # keep the last record (assuming most recent) if duplicated
    df = df[~df.index.duplicated(keep="last")]
    return df.loc[keys].reset_index().to_dict(orient="records")
