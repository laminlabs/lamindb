from typing import Union

import pandas as pd

from ...schema._table import Table
from . import exception


class SelectResult:
    """Select results."""

    def __init__(self, results: Union[list, None], model) -> None:
        if results is None:
            results = []
        self._results = results
        self._model = model

    def df(self):
        """Return list select results as a DataFrame."""
        if len(self._results) > 0:
            df = pd.DataFrame(
                [result.dict() for result in self._results],
                columns=Table.get_fields(self._model),
            )
        else:
            df = pd.DataFrame(columns=Table.get_fields(self._model))

        if "id" in df.columns:
            if "v" in df.columns:
                df = df.set_index(["id", "v"])
            else:
                df = df.set_index("id")
        return df

    def all(self):
        """Return all results of this select as a list."""
        return self._results

    def first(self):
        """Return the first result of select or None if no results are found."""
        if len(self._results) == 0:
            return None
        return self._results[0]

    def one(self):
        """Return exactly one result or raise an exception."""
        if len(self._results) == 0:
            raise exception.NoResultFound
        elif len(self._results) > 1:
            raise exception.MultipleResultsFound
        else:
            return self._results[0]

    def one_or_none(self):
        """Return at most one result or raise an exception."""
        if len(self._results) == 0:
            return None
        elif len(self._results) == 1:
            return self._results[0]
        else:
            raise exception.MultipleResultsFound
