import sqlmodel as sqm
from lndb_setup import settings


def add(*rows: sqm.SQLModel):
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
        *rows: One or multiple rows as instances of `SQLModel`.
    """
    with sqm.Session(settings.instance.db_engine()) as session:
        for row in rows:
            session.add(row)
        session.commit()
    settings.instance._update_cloud_sqlite_file()
