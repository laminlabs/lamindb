from typing import Tuple, Union

import sqlmodel as sqm
from lndb_setup import settings


def get(table: sqm.SQLModel, key: Union[str, int, Tuple]) -> Union[None, sqm.SQLModel]:
    """Retrieve an entry from a table by primary key."""
    with sqm.Session(settings.instance.db_engine()) as session:
        result = session.get(table, key)
    return result
