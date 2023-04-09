from typing import Dict, List, Tuple, Union, overload  # noqa

import sqlmodel as sqm
from lndb import settings

from .._docs import doc_args
from ._add import add, add_docs
from ._select import SelectStmt, select, select_docs


class Session:
    """Database session.

    FAQ: :doc:`/faq/session`

    It offers `.select`, `.add` and `.delete` attached to an open session, which
    is needed for lazy loading of relationships.

    The session object should be closed when it's not longer needed and
    typically used within a `with` statement::

        with Session() as ss:
            file = ss.select(ln.File, name="My test").one()
            ...
    """

    def __init__(self):
        settings.instance._cloud_sqlite_locker.lock()
        self._session = settings.instance.session()

        self._update = False

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self._session.__exit__(exc_type, exc_value, exc_tb)
        self._finish()

    def _finish(self):
        if self._update:
            settings.instance._update_cloud_sqlite_file()
        settings.instance._cloud_sqlite_locker.unlock()

    @doc_args(add_docs)
    def add(
        self, record: Union[sqm.SQLModel, List[sqm.SQLModel]], **fields
    ) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
        """{}"""  # noqa
        fields["session"] = self._session
        self._update = True
        return add(record=record, **fields)

    @doc_args(select_docs)
    def select(self, *entity: sqm.SQLModel, **fields) -> SelectStmt:
        """{}"""
        fields["session"] = self._session
        return select(*entity, **fields)

    def close(self):
        """Close the session."""
        self._session.close()
        self._finish()
