from lndb_setup import settings
from sqlmodel import Session

from .. import schema


class update:
    """Update data."""

    @classmethod
    def biometa(
        cls,
        id,
        **kwargs,
    ):
        """Query biometa by its id to update the entry fields."""
        with Session(settings.instance.db_engine()) as session:
            biometa = session.get(schema.biolab.biometa, id)
            for k, v in kwargs.items():
                biometa.__setattr__(k, int(v))
            session.add(biometa)
            session.commit()
            session.refresh(biometa)
