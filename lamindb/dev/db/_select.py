import pandas as pd
import sqlmodel as sqm
from lndb_setup import settings
from lndb_setup._settings_store import InstanceSettingsStore
from sqlalchemy.orm.exc import MultipleResultsFound, NoResultFound
from sqlmodel.main import SQLModelMetaclass


def select(*entity: sqm.SQLModel, **fields) -> "SelectStmt":
    """Select data.

    Guide: :doc:`/db/guide/select-load`.

    Returns a :class:`~lamindb.dev.db.SelectStmt` object.

    Args:
        entity: Table, tables, or tables including column specification.
        fields: Fields and values passed as keyword arguments.
    """
    # _settings_store is an internal, non-user facing variable
    _settings_store = None
    if len(fields) > 0:
        for k, v in fields.items():
            if isinstance(v, InstanceSettingsStore):
                _settings_store = fields.pop(k)
                break
    # continue with user-facing variables
    if len(entity) > 1 and len(fields) > 0:
        raise RuntimeError("Can only pass fields for a single entity.")
    elif len(fields) > 0:
        # was in `get` before, but there it leads to an inhomogeneous return type
        conditions = [getattr(entity[0], k) == v for k, v in fields.items()]
        return SelectStmt(*entity, _settings_store=_settings_store).where(*conditions)
    return SelectStmt(*entity, _settings_store=_settings_store)


def to_df(*entities, result):
    if len(result) == 0:
        if len(entities) == 1:
            entity = entities[0]
            return pd.DataFrame(columns=entity.__fields__)
        else:
            # not implemented
            return pd.DataFrame()
    row = result[0]
    if len(entities) == 1:
        entity = entities[0]
        if isinstance(entity, SQLModelMetaclass):
            records = [row.dict() for row in result]
            df = pd.DataFrame.from_records(
                records, columns=entity.__fields__, coerce_float=True
            )
            if "id" in df.columns:
                if "v" in df.columns:
                    df = df.set_index(["id", "v"])
                else:
                    df = df.set_index("id")
            return df
        elif hasattr(entity, "class_") and isinstance(entity.class_, SQLModelMetaclass):
            return pd.DataFrame(result)
        else:
            raise RuntimeError(
                "Unknown entities. Pass either a SQLModel or a field of a SQLModel."
            )
    if len(entities) > 1:
        records = []
        keys_types = []
        for entity in entities:
            if isinstance(entity, SQLModelMetaclass):
                keys_types.append((entity.__table__.name, 0))
            elif hasattr(entity, "class_") and isinstance(
                entity.class_, SQLModelMetaclass
            ):
                keys_types.append((entity.class_.__table__.name, 1))
            else:
                raise RuntimeError(
                    "Unknown entities. Pass either a SQLModel or a field of a SQLModel."
                )
        for row in result:
            record = {}
            for i, (key, type) in enumerate(keys_types):
                if type == 0:
                    sub_record = {
                        (key, field): getattr(row[i], field)
                        for field in row[i].__fields__
                    }
                    record.update(sub_record)
                else:
                    k = entities[i].expression.name
                    record[(key, k)] = row[i]
            records.append(record)
        df = pd.DataFrame(records)
        df.columns = pd.MultiIndex.from_tuples(df.columns)
        return df


class ExecStmt:
    """Executable statement."""

    def __init__(self, *, tables, stmt, _settings_store: InstanceSettingsStore = None):
        self._stmt = stmt
        self._tables = tables
        self._settings_store = _settings_store
        self._result = None

    def _execute(self):
        # cache the query result for the lifetime of the object
        if self._result is None:
            if self._settings_store is not None:
                session = sqm.Session(
                    settings._instance(self._settings_store).db_engine(),
                    expire_on_commit=False,
                )
            else:
                session = settings.instance.session()
            self._result = session.exec(self._stmt).all()
            if settings.instance._session is None:
                session.close()

    def all(self):
        """Return all result as a list."""
        self._execute()
        return self._result

    def df(self):
        """Return result as a DataFrame.

        The DataFrame will have MultiIndex columns if multiple entities are passed.
        """
        self._execute()
        return to_df(*self._tables, result=self._result)

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

    def __init__(
        self, *tables, stmt=None, _settings_store: InstanceSettingsStore = None
    ) -> None:
        self._tables = tables
        self._settings_store = _settings_store
        if stmt is None:
            stmt = sqm.select(*tables)
        super().__init__(tables=tables, stmt=stmt, _settings_store=_settings_store)

    def where(self, *conditions):
        """Pass one or multiple conditions.

        If multiple conditions are passed, they are combined using AND.

        If OR is desired, use `sqlmodel.or_`.
        """
        return SelectStmt(
            *self._tables,
            stmt=self._stmt.where(*conditions),
            _settings_store=self._settings_store
        )

    def join(self, *expression, **fields):
        """Pass a target table as an expression."""
        stmt = self._stmt.join(*expression)
        if len(fields) > 0 and len(expression) == 1:
            conditions = [getattr(expression[0], k) == v for k, v in fields.items()]
            stmt = stmt.where(*conditions)
        elif len(fields) > 0 and len(expression) > 1:
            raise RuntimeError("Can only pass fields for a single entity.")
        return SelectStmt(
            *self._tables, stmt=stmt, _settings_store=self._settings_store
        )

    def order_by(self, expression):
        """Pass a field."""
        return SelectStmt(
            *self._tables,
            stmt=self._stmt.order_by(expression),
            _settings_store=self._settings_store
        )

    def offset(self, n):
        """Pass an integer."""
        return SelectStmt(
            *self._tables,
            stmt=self._stmt.offset(n),
            _settings_store=self._settings_store
        )

    def limit(self, n):
        """Pass an integer."""
        return SelectStmt(
            *self._tables,
            stmt=self._stmt.limit(n),
            _settings_store=self._settings_store
        )
