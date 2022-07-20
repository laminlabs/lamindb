import sqlmodel as sqm
from lamindb_schema.id import id_dobject, id_gene, id_user  # noqa

import lamindb as db
from lamindb import _setup

from ..._logger import logger
from ..._setup import load_settings
from . import get_engine


class insert_if_not_exists:
    """Insert data if it does not yet exist."""

    @classmethod
    def user(cls, user_email, user_id):
        df_user = db.do.load("user")
        if user_id not in df_user.index:
            user_id = insert.user(user_email, user_id)  # type: ignore
        return user_id


class insert:
    """Insert data."""

    @classmethod
    def schema_version(cls, version, user_id):
        """User."""
        engine = get_engine()

        with sqm.Session(engine) as session:
            user = db.schema.schema_version(id=version, user_id=user_id)
            session.add(user)
            session.commit()

        load_settings()._update_cloud_sqlite_file()

    @classmethod
    def user(cls, user_email, user_id):
        """User."""
        engine = get_engine()

        with sqm.Session(engine) as session:
            user = db.schema.user(id=user_id, email=user_email)
            session.add(user)
            session.commit()
            session.refresh(user)

        load_settings()._update_cloud_sqlite_file()

        return user.id

    @classmethod
    def dobject(
        cls,
        *,
        name: str,
        file_suffix: str = None,
        interface_id: str,
        interface_v: str,
        interface_name: str,
        interface_type: str,
        dobject_id: str = None,
        dobject_v: str = "1",
    ):
        """Data object with its origin."""
        engine = get_engine()
        settings = _setup.load_settings()

        df_interface = db.do.load("interface")
        if interface_id not in df_interface.index:
            with sqm.Session(engine) as session:
                interface = db.schema.interface(
                    id=interface_id,
                    v=interface_v,
                    name=interface_name,
                    type=interface_type,
                    user_id=settings.user_id,
                )
                session.add(interface)
                session.commit()
            logger.info(
                f"Added notebook {interface_name!r} ({interface_id}, {interface_v}) by"
                f" user {settings.user_email} ({settings.user_id})."
            )

        with sqm.Session(engine) as session:
            dobject = db.schema.dobject(
                id=dobject_id,
                v=dobject_v,
                name=name,
                interface_id=interface_id,
                interface_v=interface_v,
                file_suffix=file_suffix,
            )
            session.add(dobject)
            session.commit()
            session.refresh(dobject)

        load_settings()._update_cloud_sqlite_file()

        return dobject.id

    @classmethod
    def gene(
        cls,
        name: str,
        feature_type: str = None,
        database: str = None,
        species: str = None,
        curated: bool = None,
    ):
        """Genes in a dobject."""
        engine = get_engine()

        with sqm.Session(engine) as session:
            gene = db.schema.gene(
                id=id_gene(),
                name=name,
                feature_type=feature_type,
                database=database,
                species=species,
                curated=curated,
            )
            session.add(gene)
            session.commit()
            session.refresh(gene)

        load_settings()._update_cloud_sqlite_file()

        return gene.id
