from typing import List, Union, overload  # noqa

import sqlmodel as sqm
from lndb_setup import settings


@overload
def add(record: sqm.SQLModel) -> sqm.SQLModel:
    ...


# Currently seeing the following error without type ignore:
# Overloaded function signature 2 will never be matched: signature 1's parameter
# type(s) are the same or broader
@overload
def add(records: List[sqm.SQLModel]) -> List[sqm.SQLModel]:  # type: ignore
    ...


def add(  # type: ignore
    record: Union[sqm.SQLModel, List[sqm.SQLModel]]
) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
    """Insert or update records.

    Inserts a new record if the corresponding row doesn't exist.
    Updates the corresponding row with the record if it exists.

    To update a row, query it with `.get` or `.select` and modify it before
    passing it to `add`.

    Guide: :doc:`/db/guide/add-delete`.

    Example:

    >>> # insert a new record
    >>> db.add(wetlab.experiment(name="My test", biometa_id=test_id))
    >>> # update an existing record
    >>> experiment = ln.get(wetlab.experiment, experiment_id)
    >>> experiment.name = "New name"
    >>> db.add(experiment)

    Args:
        record: One or multiple records as instances of `SQLModel`.
    """
    if isinstance(record, list):
        records = record
    else:
        records = [record]
    with sqm.Session(settings.instance.db_engine()) as session:
        for record in records:
            session.add(record)
        session.commit()
        for record in records:
            session.refresh(record)
    settings.instance._update_cloud_sqlite_file()
    if len(records) > 1:
        return records
    else:
        return records[0]
