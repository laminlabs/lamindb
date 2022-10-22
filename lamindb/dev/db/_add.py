from typing import List, Union, overload  # noqa

import sqlmodel as sqm
from lndb_setup import settings


@overload
def add(rows: sqm.SQLModel) -> sqm.SQLModel:
    ...


# Currently seeing the following error without type ignore:
# Overloaded function signature 2 will never be matched: signature 1's parameter
# type(s) are the same or broader
@overload
def add(rows: List[sqm.SQLModel]) -> List[sqm.SQLModel]:  # type: ignore
    ...


def add(
    rows: Union[sqm.SQLModel, List[sqm.SQLModel]]
) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    """Insert or update rows.

    Inserts a new row if it doesn't exist. Updates the row if it exists already.

    To achieve the latter, query the `row` with `.get` or `.select` before
    passing it to add.

    Guide: :doc:`/db/guide/add-delete`.

    Example:

    >>> experiment = ln.get(wetlab.experiment, "93jIJFla")
    >>> experiment.name = "New name"
    >>> add(experiment)

    Args:
        rows: One or multiple rows as instances of `SQLModel`.
    """
    if not isinstance(rows, list):
        rows = [rows]
    with sqm.Session(settings.instance.db_engine()) as session:
        for row in rows:
            session.add(row)
        session.commit()
        for row in rows:
            session.refresh(row)
    settings.instance._update_cloud_sqlite_file()
    if len(rows) > 1:
        return rows
    else:
        return rows[0]
