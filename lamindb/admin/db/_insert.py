import sqlmodel as sqm

import lamindb as db
from lamindb import setup

from ..._logger import logger
from ...dev.id import id_dobject, id_user  # noqa
from . import get_engine


class insert_if_not_exists:
    """Insert data if it does not yet exist."""

    @classmethod
    def user(cls, user_name):
        df_user = db.do.load("user")

        if user_name in df_user.name.values:
            user_id = df_user.index[df_user.name == user_name][0]
            logger.info(f"Logged in {user_name} ({user_id}).")
        else:
            user_id = insert.user(user_name)  # type: ignore
            logger.info(f"Signed up user {user_name} ({user_id}).")

        return user_id


class insert:
    """Insert data."""

    @classmethod
    def user(cls, user_name):
        """User."""
        engine = get_engine()

        with sqm.Session(engine) as session:
            user = db.model.user(name=user_name)
            session.add(user)
            session.commit()
            session.refresh(user)

        return user.id

    @classmethod
    def dobject(
        cls, name: str, *, interface_id: str = None, interface_name: str = None
    ):
        """Data object with its origin."""
        engine = get_engine()
        settings = setup.settings()

        if interface_id is None:
            from nbproject import meta

            interface_id = meta.store.id
            interface_name = meta.live.title
            dependency_string = " ".join(
                [pkg + f"=={ver}" for pkg, ver in meta.live.dependency.items()]
            )
            interface_dependency = dependency_string
            interface_type = "nbproject"

            if interface_name is None:
                raise RuntimeError(
                    "Can only ingest from notebook with title. Please set a title!"
                )
        else:
            interface_dependency = None
            interface_type = "other"

        df_interface = db.do.load("interface")
        if interface_id not in df_interface.index:
            with sqm.Session(engine) as session:
                # can remove underscore once we explicitly
                # migrate to _id suffixes for id columns
                interface = db.model.interface(
                    id=interface_id,
                    name=interface_name,
                    dependency=interface_dependency,
                    type=interface_type,
                    user_id=settings.user_id,
                )
                session.add(interface)
                session.commit()
            logger.info(
                f"Added notebook {interface_name!r} ({interface_id}) by user"
                f" {settings.user_name} ({settings.user_id})."
            )

        with sqm.Session(engine) as session:
            dobject = db.model.dobject(name=name, interface_id=interface_id)
            session.add(dobject)
            session.commit()
            session.refresh(dobject)

        return dobject.id
