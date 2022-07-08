import sqlmodel as sqm

import lamindb as db
from lamindb import _setup

from ..._logger import logger
from ..._setup import load_settings
from ...dev.id import id_dobject, id_user  # noqa
from . import get_engine


class insert_if_not_exists:
    """Insert data if it does not yet exist."""

    @classmethod
    def user(cls, user_email, user_id):
        df_user = db.do.load("user")

        if user_id in df_user.index:
            logger.info(f"Instance DB: Logged in {user_email} ({user_id}).")
        else:
            user_id = insert.user(user_email, user_id)  # type: ignore
            logger.info(f"Instance DB: Signed up {user_email} ({user_id}).")

        return user_id


class insert:
    """Insert data."""

    @classmethod
    def user(cls, user_email, user_id):
        """User."""
        engine = get_engine()

        with sqm.Session(engine) as session:
            user = db.model.user(id=user_id, email=user_email)
            session.add(user)
            session.commit()
            session.refresh(user)

        load_settings()._update_cloud_sqlite_file()

        return user.id

    @classmethod
    def dobject(
        cls,
        name: str,
        suffix: str = None,
        *,
        interface_id: str = None,
        interface_name: str = None,
    ):
        """Data object with its origin."""
        engine = get_engine()
        settings = _setup.load_settings()

        if interface_id is None:
            from nbproject import meta

            interface_id = meta.store.id
            interface_name = meta.live.title
            interface_type = "nbproject"

            if interface_name is None:
                raise RuntimeError(
                    "Can only ingest from notebook with title. Please set a title!"
                )
        else:
            interface_type = "other"

        df_interface = db.do.load("interface")
        if interface_id not in df_interface.index:
            with sqm.Session(engine) as session:
                # can remove underscore once we explicitly
                # migrate to _id suffixes for id columns
                interface = db.model.interface(
                    id=interface_id,
                    name=interface_name,
                    type=interface_type,
                    user_id=settings.user_id,
                )
                session.add(interface)
                session.commit()
            logger.info(
                f"Added notebook {interface_name!r} ({interface_id}) by user"
                f" {settings.user_email} ({settings.user_id})."
            )

        with sqm.Session(engine) as session:
            dobject = db.model.dobject(
                name=name, interface_id=interface_id, suffix=suffix
            )
            session.add(dobject)
            session.commit()
            session.refresh(dobject)

        load_settings()._update_cloud_sqlite_file()

        return dobject.id
