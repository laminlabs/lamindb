from typing import Union

from lndb_setup import settings
from sqlmodel import Session, select
from sqlmodel.sql.expression import Select, SelectOfScalar

from .. import schema


class query:
    """Query literal (semantic) data."""

    @classmethod
    def id(cls, entity_name: str, id: Union[str, tuple]):
        """Query a single row by its id column with the primary key."""
        engine = settings.instance.db_engine()
        with Session(engine) as session:
            for module in [schema.core, schema.biolab, schema.bionty]:
                try:
                    return session.get(getattr(module, entity_name), id)
                except AttributeError:
                    continue

    @classmethod
    def readout_type(cls, name: str = None, platform: str = None):
        """Query from the readout_type table."""
        with Session(settings.instance.db_engine()) as session:
            stmt = select(schema.biolab.readout_type).where(
                schema.biolab.readout_type.name == name, platform == platform
            )
            # Will remove after this is fixed:
            # https://github.com/tiangolo/sqlmodel/pull/234
            SelectOfScalar.inherit_cache = True  # type: ignore
            Select.inherit_cache = True  # type: ignore
            results = session.exec(stmt).all()

        return results

    @classmethod
    def species(cls, common_name: str, taxon_id: str = None):
        """Query from the species table."""
        # if taxon_id is provided, query by taxon_id
        if taxon_id is not None:
            stmt = select(schema.bionty.species).where(
                schema.bionty.species.taxon_id == taxon_id
            )
        else:
            stmt = select(schema.bionty.species).where(
                schema.bionty.species.common_name == common_name
            )
        with Session(settings.instance.db_engine()) as session:
            # Will remove after this is fixed:
            # https://github.com/tiangolo/sqlmodel/pull/234
            SelectOfScalar.inherit_cache = True  # type: ignore
            Select.inherit_cache = True  # type: ignore
            results = session.exec(stmt).one()

        return results

    @classmethod
    def dobject_biometa(cls, dobject_id: str):
        """Query from the readout_type table."""
        with Session(settings.instance.db_engine()) as session:
            stmt = select(schema.core.dobject_biometa).where(
                schema.core.dobject_biometa.dobject_id == dobject_id
            )
            # Will remove after this is fixed:
            # https://github.com/tiangolo/sqlmodel/pull/234
            SelectOfScalar.inherit_cache = True  # type: ignore
            Select.inherit_cache = True  # type: ignore
            results = session.exec(stmt).all()

        return results
