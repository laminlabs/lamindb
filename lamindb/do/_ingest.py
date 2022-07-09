from pathlib import Path
from typing import Union

import sqlmodel as sqm
from loguru import logger

import lamindb as db
from lamindb._setup import load_settings

from ..admin.db import get_engine
from ..dev.file import store_file


def track_ingest(dobject_id):
    engine = get_engine()

    from nbproject import meta

    settings = load_settings()

    user_id = settings.user_id

    interface_id = meta.store.id

    with sqm.Session(engine) as session:
        track_do = db.model.track_do(
            type="ingest",
            user_id=user_id,
            interface_id=interface_id,
            dobject_id=dobject_id,
        )
        session.add(track_do)
        session.commit()
        session.refresh(track_do)

    settings._update_cloud_sqlite_file()

    return track_do.id


def ingest(filepath, integrity: Union[bool, None] = None):
    """Ingest file.

    Args:
        filepath: The filepath.
        integrity: Check the integrity of the notebook.

    We primarily work with base62 IDs.

    ====== =========
    len_id n_entries
    ====== =========
    1      >6e+01
    2      >4e+03
    3      >2e+05
    4      >1e+07
    5      >9e+08
    6      >6e+10
    7      >4e+12
    8      >2e+14
    9      >1e+16
    12     >3e+21 (nbproject id)
    20     >7e+35 (~UUID)
    ====== =========
    """
    from nbproject import meta, publish

    settings = load_settings()
    storage_dir = settings.storage_dir

    storage_dir = Path(storage_dir)

    filepath = Path(filepath)

    from lamindb.admin.db import insert

    dobject_id = insert.dobject(filepath.stem, filepath.suffix)

    dobject_storage_key = f"{dobject_id}{filepath.suffix}"
    store_file(filepath, dobject_storage_key)

    track_ingest(dobject_id)

    logger.info(
        f"Added file {filepath.name} ({dobject_id}) from notebook"
        f" {meta.live.title!r} ({meta.store.id}) by user"
        f" {settings.user_email} ({settings.user_id})."
    )

    from nbproject import meta

    if integrity is None:
        if meta._env == "lab":
            integrity = True
        else:
            integrity = False
            logger.warning(
                "Consider using Jupyter Lab for ingesting data!\n"
                "Interactive notebook integrity checks are currently only supported on Jupyter Lab.\n"  # noqa
                "Alternatively, manually save your notebook directly before calling `do.ingest(..., integrity=True)`."  # noqa
            )
    publish(integrity=integrity)
