import shutil
from pathlib import Path

import sqlmodel as sqm

import lamindb as db

from ..admin.db import get_engine


def track_ingest(file_id):

    engine = get_engine()

    from nbproject import meta

    from lamindb._configuration import user_id

    interface_id = meta.id

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
    from lamindb._configuration import storage_root

    storage_root = Path(storage_root)

    filepath = Path(filepath)
    filename = filepath.name

    from lamindb.admin.db import insert

    file_id = insert.file(filename)

    storage_name = str(filepath.stem) + f"--lndb-{file_id}" + str(filepath.suffix)

    storage_path = storage_root / storage_name

    shutil.copyfile(filepath, storage_path)

    import nbproject

    from lamindb._configuration import user_id

    track_ingest(file_id)
    print(f"added file {file_id} from source {nbproject.meta.id} by user {user_id}")
