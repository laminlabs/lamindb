import lnschema_core as core
from lamin_logger import logger
from lndb_setup import settings

from .dev._core import filepath_from_dobject
from .dev.db._track_usage import track_usage
from .dev.file import load_to_memory


def populate_runin(dobject: core.DObject, run: core.Run):
    session = settings.instance.session()
    result = session.get(core.link.RunIn, (run.id, dobject.id))
    if result is None:
        session.add(
            core.link.RunIn(
                run_id=run.id,
                dobject_id=dobject.id,
            )
        )
        session.commit()
        logger.info(f"Added dobject ({dobject.id}) as input for run ({run.id}).")
        settings.instance._update_cloud_sqlite_file()
    if settings.instance._session is None:
        session.close()


def load(dobject: core.DObject, stream: bool = False):
    """Load data object.

    Returns object associated with the stored `dobject`.

    Populates `RunIn` when called from a notebook.

    Guide: :doc:`/db/guide/select-load`.
    """
    if stream and dobject.suffix not in (".h5ad", ".zarr"):
        logger.warning(f"Ignoring stream option for a {dobject.suffix} object.")

    filepath = filepath_from_dobject(dobject)
    from lamindb._nb import _run as nb_run

    if nb_run is None:
        logger.warning(
            "Input tracking for runs through `load` is currently only implemented for"
            " notebooks."
        )
    else:
        populate_runin(dobject, nb_run)
    track_usage(dobject.id, "load")
    return load_to_memory(filepath, stream=stream)
