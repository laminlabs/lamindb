import sqlmodel as sqm
from lamindb_schema.id import id_dobject, id_user  # noqa

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
        jupynb_id: str,
        jupynb_v: str,
        jupynb_name: str,
        jupynb_type: str,
        dobject_id: str = None,
        dobject_v: str = "1",
    ):
        """Data object with its origin."""
        engine = get_engine()
        settings = _setup.load_settings()

        df_jupynb = db.do.load("jupynb")
        if jupynb_id not in df_jupynb.index:
            with sqm.Session(engine) as session:
                jupynb = db.schema.jupynb(
                    id=jupynb_id,
                    v=jupynb_v,
                    name=jupynb_name,
                    type=jupynb_type,
                    user_id=settings.user_id,
                )
                session.add(jupynb)
                session.commit()
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {settings.user_email} ({settings.user_id})."
            )

        with sqm.Session(engine) as session:
            dobject = db.schema.dobject(
                id=dobject_id,
                v=dobject_v,
                name=name,
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
                file_suffix=file_suffix,
            )
            session.add(dobject)
            session.commit()
            session.refresh(dobject)

        load_settings()._update_cloud_sqlite_file()

        return dobject.id
