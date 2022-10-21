import pandas as pd
import sqlmodel as sqm
from lndb_setup import settings
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound


def select(*tables_or_columns: sqm.SQLModel):
    """Select rows.

    Guide: :doc:`/db/guide/select-load`.

    Returns a :class:`~lamindb.dev.db.SelectStmt` object.

    Args:
       tables: Tables or columns.
    """
    return SelectStmt(*tables_or_columns)


class ExecStmt:
    """Executable statement."""

    def __init__(self, *, tables, stmt):
        self._stmt = stmt
        self._tables = tables
        self._result = None

    def _execute(self):
        # cache the query result for the lifetime of the object
        if self._result is None:
            with sqm.Session(settings.instance.db_engine()) as session:
                self._result = session.exec(self._stmt).all()

    def all(self):
        """Return all result as a list."""
        self._execute()
        return self._result

    def df(self):
        """Return result as a DataFrame or list of DataFrames."""
        self._execute()
        dfs = []
        for itable, table in enumerate(self._tables):
            if len(self._result) > 0:
                df = pd.DataFrame(
                    [
                        result.dict()
                        if len(self._tables) == 1
                        else result[itable].dict()
                        for result in self._result
                    ],
                    columns=table.__fields__,
                )
            else:
                df = pd.DataFrame(columns=table.__fields__)
            dfs.append(df)

        # turn id and v columns into the index
        for i, df in enumerate(dfs):
            if "id" in df.columns:
                if "v" in df.columns:
                    dfs[i] = df.set_index(["id", "v"])
                else:
                    dfs[i] = df.set_index("id")

        if len(dfs) == 1:
            return dfs[0]
        else:
            return dfs

    def first(self):
        """Return the first result of select or None if no result are found."""
        self._execute()
        if len(self._result) == 0:
            return None
        return self._result[0]

    def one(self):
        """Return exactly one result or raise an exception."""
        self._execute()
        if len(self._result) == 0:
            raise NoResultFound
        elif len(self._result) > 1:
            raise MultipleResultsFound
        else:
            return self._result[0]

    def one_or_none(self):
        """Return at most one result or raise an exception."""
        self._execute()
        if len(self._result) == 0:
            return None
        elif len(self._result) == 1:
            return self._result[0]
        else:
            raise MultipleResultsFound


class SelectStmt(ExecStmt):
    """Executable select statement.

    Offers an exact subset of SQLAlchemy of constraining select with `where`,
    `join`, `order_by`, `offset`, `limit`.

    Offers accessing the results of the executed statement via `all`, `df`,
    `one`, `one_or_none` through the base class `ExecStmt`.
    """

    def __init__(self, *tables, stmt=None) -> None:
        self._tables = tables
        if stmt is None:
            stmt = sqm.select(*tables)
        super().__init__(tables=tables, stmt=stmt)

    def where(self, *conditions):
        """Pass one or multiple conditions.

        If multiple conditions are passed, they are combined using AND.

        If OR is desired, use `sqlmodel.or_`.
        """
        return SelectStmt(*self._tables, stmt=self._stmt.where(*conditions))

    def join(self, *expression, **kwargs):
        """Pass one or multiple conditions.

        If multiple conditions are passed, they are combined using AND.

        If OR is desired, use `sqlmodel.or_`.
        """
        return SelectStmt(*self._tables, stmt=self._stmt.join(*expression, **kwargs))

    def order_by(self, expression):
        """Pass one or multiple conditions.

        If multiple conditions are passed, they are combined using AND.

        If OR is desired, use `sqlmodel.or_`.
        """
        return SelectStmt(*self._tables, stmt=self._stmt.order_by(expression))

    def offset(self, n):
        """Pass an integer."""
        return SelectStmt(*self._tables, stmt=self._stmt.offset(n))

    def limit(self, n):
        """Pass an integer."""
        return SelectStmt(*self._tables, stmt=self._stmt.limit(n))
