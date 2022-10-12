import sqlmodel as sqm
from lndb_setup import settings


def session() -> sqm.Session:
    """Get connection session to DB engine.

    Returns a `sqlmodel.Session` object.
    """
    return sqm.Session(settings.instance.db_engine())
