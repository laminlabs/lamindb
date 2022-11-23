import os

import sqlmodel as sqm
from lndb_setup import settings

_session = None
if "LAMIN_SKIP_MIGRATION" not in os.environ:
    _session = sqm.Session(settings.instance.db_engine(), expire_on_commit=False)


def session() -> sqm.Session:
    """Get connection session to DB engine.

    Returns a `sqlmodel.Session` object.
    """
    if _session is not None:
        assert _session
        return _session
    else:
        return sqm.Session(settings.instance.db_engine(), expire_on_commit=False)
