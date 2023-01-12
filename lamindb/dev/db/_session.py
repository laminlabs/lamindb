from typing import Dict, List, Tuple, Union, overload  # noqa

import sqlmodel as sqm
from lndb_setup import settings

from .._docs import doc_args
from ._add import add, add_docs


class Session:
    """Database session.

    A class that allows to work with records that have an open connection to the
    database.

    This is needed for lazy loading of relationships.
    """

    def __init__(self):
        self._session = settings.instance.session()

    @doc_args(add_docs)
    def add(
        self,
        record: Union[sqm.SQLModel, List[sqm.SQLModel]],
        use_fsspec: bool = True,
        **fields
    ) -> Union[sqm.SQLModel, List[sqm.SQLModel]]:
        """{add_docs}"""  # noqa
        fields["session"] = self._session
        return add(record=record, use_fsspec=use_fsspec, **fields)

    def close(self):
        """Close the session."""
        self._session.close()
