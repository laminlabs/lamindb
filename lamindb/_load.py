import lnschema_core as core
import sqlmodel as sqm
from lamin_logger import logger
from lndb_setup import settings
from nbproject import meta

from .dev._core import filepath_from_dobject
from .dev.db._track_usage import track_usage
from .dev.file import load_to_memory


def populate_run_in(dobject):
    jupynb_id = meta.store.id
    jupynb_v = meta.store.version  # version to be set in publish()
    jupynb_name = meta.live.title
    engine = settings.instance.db_engine()

    committed = False

    with sqm.Session(engine) as session:
        result = session.get(core.Jupynb, (jupynb_id, jupynb_v))
        if result is None:
            session.add(
                core.Jupynb(
                    id=jupynb_id, v=jupynb_v, name=jupynb_name, user_id=settings.user.id
                )
            )
            session.commit()
            run = core.Run(
                jupynb_id=jupynb_id,
                jupynb_v=jupynb_v,
            )
            session.add(run)
            session.commit()
            committed = True
            logger.info(
                f"Added notebook {jupynb_name!r} ({jupynb_id}, {jupynb_v}) by"
                f" user {settings.user.handle}."
            )
        else:
            run = session.exec(
                sqm.select(core.Run).where(
                    core.Run.jupynb_id == jupynb_id,
                    core.Run.jupynb_v == jupynb_v,
                )
            ).one()
        result = session.get(core.RunIn, (run.id, dobject.id))
        if result is None:
            session.add(
                core.RunIn(
                    run_id=run.id,
                    dobject_id=dobject.id,
                )
            )
            session.commit()
            committed = True
            logger.info(f"Added dobject ({dobject.id}) as input for run ({run.id}).")
    if committed:
        # nothing to update if the db file wasn't changed
        settings.instance._update_cloud_sqlite_file()


def load(dobject: core.DObject, stream: bool = False):
    """Load data object.

    Returns object associated with the stored `dobject`.

    Populates `run_in`.

    Guide: :doc:`/db/guide/select-load`.
    """
    if stream and dobject.suffix not in (".h5ad", ".zarr"):
        logger.warning(f"Ignoring stream option for a {dobject.suffix} object.")

    filepath = filepath_from_dobject(dobject)
    populate_run_in(dobject)
    track_usage(dobject.id, "load")
    return load_to_memory(filepath, stream=stream)
