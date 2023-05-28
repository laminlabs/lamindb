from typing import Dict, List, TypeVar, Union

import numpy as np
import pandas as pd
from lamin_logger import logger
from sqlalchemy.orm.attributes import InstrumentedAttribute
from sqlmodel import SQLModel

from .dev.db._select import select

ListLike = TypeVar("ListLike", pd.Series, list, np.array)


def parse(
    iterable: Union[ListLike, pd.DataFrame],
    field: Union[InstrumentedAttribute, Dict[str, InstrumentedAttribute]],
    from_bionty: bool = True,
) -> List[SQLModel]:
    """Parse a dataset column based on a SQLModel entity field.

    Guide: :doc:`/guide/ontologies`.

    Args:
        iterable: a `ListLike` of values or a `DataFrame`.
        field: if iterable is `ListLike`: a `SQLModel` field to parse into.
               if iterable is `DataFrame`: a dict of {column_name : SQLModel_field}.
        from_bionty: whether to auto-complete the fields of bionty tables.
                     only effective if iterable is `ListLike`.

    Returns:
        A list of SQLModel records.
    """
    if isinstance(iterable, pd.DataFrame):
        if not isinstance(field, dict):
            raise TypeError(
                "field must be a dictionary of {column_name : SQLModel_field}!"
            )
        column_mapper = {colname: f.name for colname, f in field.items()}
        class_mapper = {f.name: f.class_ for f in field.values()}
        if len(set(class_mapper.values())) > 1:
            raise NotImplementedError("fields must from the same entity!")
        entity = list(class_mapper.values())[0]
        df = iterable.rename(columns=column_mapper)
        df = df.dropna()
        iterable = df.to_dict(orient="records")
    else:
        if not isinstance(field, InstrumentedAttribute):
            raise TypeError("field must be a SQLModel field!")
        entity = field.class_
        iterable = set(iterable)

    records = []
    for i in iterable:
        # No entries are made for NAs, '', None
        if pd.isnull(i) or i == "":
            continue
        if isinstance(i, dict):
            kwargs = i
        else:
            kwargs = {field.name: i}  # type:ignore
        record = _create_record(entity, **kwargs, from_bionty=from_bionty)
        records.append(record)
    return records


def _create_record(entity: SQLModel, from_bionty: bool = True, **kwargs):
    """Create a record from bionty with configured ontologies."""
    if len(kwargs) > 1:
        from_bionty = False
    kwargs = {k: v for k, v in kwargs.items() if k in entity.__fields__}
    record = select(entity, **kwargs).one_or_none()
    if record is None:
        if from_bionty:
            try:
                # bionty tables
                record = entity.from_bionty(**kwargs)
            except (ValueError, KeyError):
                logger.warning(
                    f"No entry found in bionty with {kwargs}!"
                    " Couldn't populate additional fields..."
                )
                record = entity(**kwargs)
            except AttributeError:
                # no-bionty tables
                record = entity(**kwargs)
        else:
            record = entity(**kwargs)
    return record
