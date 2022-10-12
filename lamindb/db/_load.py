import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings
from nbproject import meta

from ..dev._core import filepath_from_dobject
from ..dev.db._track_usage import track_usage
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
            session.commit()
            dtransform = core.dtransform(
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
            )
            session.add(dtransform)
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
        result = session.get(core.dtransform_in, (dtransform.id, dobject.id))
        if result is None:
            session.add(
                core.dtransform_in(
                    dtransform_id=dtransform.id,
                    dobject_id=dobject.id,
                )
            )
            session.commit()
            committed = True
            logger.info(
                f"Added dobject ({dobject.id}) as input for dtransform"
                f" ({dtransform.id})."
            )
    if committed:
        # nothing to update if the db file wasn't changed
        settings.instance._update_cloud_sqlite_file()


def load(dobject: core.dobject, stream: bool = False):
    """Load data object into memory.

    Returns object associated with the stored `dobject`.

    Populates `dtransform_in`.

    Guide: :doc:`/db/guide/select-load`.
    """
    if stream and dobject.suffix not in (".h5ad", ".zarr"):
        logger.warning(f"Ignoring stream option for a {dobject.suffix} object.")

    filepath = filepath_from_dobject(dobject)
    populate_dtransform_in(dobject)
    track_usage(dobject.id, "load")
    return load_to_memory(filepath, stream=stream)
