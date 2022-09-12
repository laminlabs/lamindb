import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings
from lnschema_core import id
from nbproject import meta

from ..dev import filepath_from_dobject, track_usage
from ..dev.file import load_to_memory
from ..schema import core


def populate_dtransform_in(dobject):
    jupynb_id = meta.store.id
    jupynb_v = meta.store.version  # version to be set in publish()
    jupynb_name = meta.live.title
    engine = settings.instance.db_engine()

    committed = False

    with sqm.Session(engine) as session:
        result = session.get(core.jupynb, (jupynb_id, jupynb_v))
        if result is None:
            session.add(
                core.jupynb(
                    id=jupynb_id, v=jupynb_v, name=jupynb_name, user_id=settings.user.id
                )
            )
            dtransform_id = id.id_dtransform()
            session.add(
                core.dtransform(
                    id=dtransform_id,
                    jupynb_id=jupynb_id,
                    jupynb_v=jupynb_v,
                )
            )
            session.commit()
            committed = True
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {settings.user.handle}."
            )
        else:
            dtransform = session.exec(
                sqm.select(core.dtransform).where(
                    core.dtransform.jupynb_id == jupynb_id,
                    core.dtransform.jupynb_v == jupynb_v,
                )
            ).one()
            dtransform_id = dtransform.id
        result = session.get(core.dtransform_in, (dtransform_id, dobject.id, dobject.v))
        if result is None:
            session.add(
                core.dtransform_in(
                    dtransform_id=dtransform_id,
                    dobject_id=dobject.id,
                    dobject_v=dobject.v,
                )
            )
            session.commit()
            committed = True
            logger.info(
                f"Added dobject ({dobject.id}, {dobject.v}) as input for dtransform"
                f" ({dtransform_id})."
            )
    if committed:
        # nothing to update if the db file wasn't changed
        settings.instance._update_cloud_sqlite_file()


def load(dobject: core.dobject):
    """Load `dobject` into memory.

    Returns object associated with the stored `dobject`.

    Populates `dtransform_in`.
    """
    filepath = filepath_from_dobject(dobject)
    populate_dtransform_in(dobject)
    track_usage(dobject.id, dobject.v, "load")
    return load_to_memory(filepath)
