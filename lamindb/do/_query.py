from typing import Union

import pandas as pd
from sqlmodel import Session, select
from sqlmodel.sql.expression import Select, SelectOfScalar

from .. import schema
from ..dev.db import get_engine


class query:
    """Query literal (semantic) data."""

    @classmethod
    def id(cls, entity_name: str, id: Union[str, tuple]):
        """Query a single row by its id column with the primary key."""
        engine = get_engine()
        with Session(engine) as session:
            return session.get(getattr(schema.core, entity_name), id)

    @classmethod
    def readout_type(cls, name: str = None, platform: str = None):
        """Query from the readout_type table."""
        engine = get_engine(future=False)
        statement = select(schema.biolab.readout_type).where(
            schema.biolab.readout_type.name == name, platform == platform
        )

        # Will remove after this is fixed:
        # https://github.com/tiangolo/sqlmodel/pull/234
        SelectOfScalar.inherit_cache = True  # type: ignore
        Select.inherit_cache = True  # type: ignore

        results = pd.read_sql_query(statement, engine)

        return results
