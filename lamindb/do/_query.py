from typing import Union

from sqlmodel import Session

from .. import schema
from ..dev.db import get_engine


class query:
    """Query literal (semantic) data."""

    @classmethod
    def id(cls, entity_name: str, id: Union[str, tuple]):
        """Query a single row by its id column with the primary key."""
        engine = get_engine()
        with Session(engine) as session:
            return session.get(getattr(schema.provenance, entity_name), id)
