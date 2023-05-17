from typing import Callable, List, TypeVar

import numpy as np
import pandas as pd
from lamin_logger import logger
from sqlalchemy.orm.attributes import InstrumentedAttribute
from sqlmodel import SQLModel

from .dev.db._select import select

ListLike = TypeVar("ListLike", pd.Series, list, np.array)
SQLModelField = TypeVar("SQLModelField", Callable, InstrumentedAttribute)


def parse(
    iterable: ListLike, field: SQLModelField, from_bionty: bool = True
) -> List[SQLModel]:
    """Parse a dataset column based on a SQLModel entity field.

    Guide: :doc:`/guide/ontologies`.

    Args:
        iterable: a `ListLike` of values.
        field: a `SQLModel` field that correspond to the input list.
        from_bionty: whether to use bionty for auto-completion of the fields.

    Returns:
        A list of SQLModel records.
    """
    records = []
    entity = field.class_  # type:ignore
    for i in set(iterable):
        # No entries are made for NAs, '', None
        if pd.isnull(i) or i == "":
            continue
        kwargs = {field.name: i}  # type:ignore
        record = _create_record(entity, **kwargs, from_bionty=from_bionty)
        records.append(record)
    return records


def _create_record(entity: SQLModel, from_bionty: bool = True, **kwargs):
    """Create a record from bionty with configured ontologies."""
    record = select(entity, **kwargs).one_or_none()
    if record is None:
        if from_bionty:
            try:
                record = entity.from_bionty(**kwargs)
            except AttributeError:
                record = entity(**kwargs)
            except (ValueError, KeyError):
                logger.warning(
                    f"No entry found in bionty with {kwargs}!"
                    " Couldn't populate additional fields..."
                )
                record = entity(**kwargs)
        else:
            record = entity(**kwargs)
    return record
