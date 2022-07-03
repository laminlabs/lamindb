from pathlib import Path

import sqlmodel as sqm
from loguru import logger

import lamindb as db
from lamindb import setup

from ..admin.db import get_engine
from ..dev.file import store_file


def track_ingest(file_id):
    engine = get_engine()

    from nbproject import meta

    user_id = setup.settings().user_id

    interface_id = meta.store.id

    with sqm.Session(engine) as session:
        track_do = db.model.track_do(
            type="ingest",
            user_id=user_id,
            interface_id=interface_id,
            file_id=file_id,
        )
        session.add(track_do)
        session.commit()
        session.refresh(track_do)

    return track_do.id


def ingest(filepath):
    """Ingest file.

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

    settings = setup.settings()
    storage_dir = settings.storage_dir

    storage_dir = Path(storage_dir)

    filepath = Path(filepath)
    filename = filepath.name

    from lamindb.admin.db import insert

    file_id = insert.file(filename)

    filekey = f"{file_id}.{filepath.suffix}"
    store_file(filepath, filekey)

    track_ingest(file_id)
    logger.info(
        f"Added file {file_id} from notebook {meta.live.title!r} ({meta.store.id}) by"
        f" user {settings.user_name} ({settings.user_id}).",
        flush=True,
    )
    publish(integrity=False)  # noqa
