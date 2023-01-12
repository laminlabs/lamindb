from typing import Dict, List, Tuple, Union, overload  # noqa

import sqlmodel as sqm
from lndb_setup import settings

from .._docs import doc_args
from ._add import add, add_docs
from ._select import SelectStmt, select, select_docs


class Session:
    """Database session.

    It offers `.select`, `.add` and `.delete` attached to an open session, which
    is needed for lazy loading of relationships.

    The session object should be closed when it's not longer needed and
    typically used within a `with` statement::

        with Session() as ss:
            dobject = ss.select(lns.DObject, name="My test").one()
            ...
    """

    def __init__(self):
        self._session = settings.instance.session()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_value, exc_tb):
        self._session.__exit__(exc_type, exc_value, exc_tb)

    @doc_args(add_docs)
    def add(
        self,
        record: Union[sqm.SQLModel, List[sqm.SQLModel]],
        use_fsspec: bool = True,
        **fields
    ) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
        """{}"""  # noqa
        fields["session"] = self._session
        return add(record=record, use_fsspec=use_fsspec, **fields)

    @doc_args(select_docs)
    def select(self, *entity: sqm.SQLModel, **fields) -> SelectStmt:
        """{}"""
        fields["session"] = self._session
        return select(*entity, **fields)

    def close(self):
        """Close the session."""
        self._session.close()
