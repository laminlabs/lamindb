from typing import Union

from lndb_setup import load_or_create_instance_settings
from sqlmodel import Session, select
from sqlmodel.sql.expression import Select, SelectOfScalar

from .. import schema


class query:
    """Query literal (semantic) data."""

    @classmethod
    def id(cls, entity_name: str, id: Union[str, tuple]):
        """Query a single row by its id column with the primary key."""
        settings = load_or_create_instance_settings()
        engine = settings.db_engine()
        with Session(engine) as session:
            return session.get(getattr(schema.core, entity_name), id)

    @classmethod
    def readout_type(cls, name: str = None, platform: str = None):
        """Query from the readout_type table."""
        settings = load_or_create_instance_settings()

        with Session(settings.db_engine()) as session:
            stmt = select(schema.biolab.readout_type).where(
                schema.biolab.readout_type.name == name, platform == platform
            )
            # Will remove after this is fixed:
            # https://github.com/tiangolo/sqlmodel/pull/234
            SelectOfScalar.inherit_cache = True  # type: ignore
            Select.inherit_cache = True  # type: ignore
            results = session.exec(stmt).all()

        return results
