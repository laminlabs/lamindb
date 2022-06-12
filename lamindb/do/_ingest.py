import shutil
from pathlib import Path

import sqlalchemy as sql

from ..admin.db import get_engine
from ..dev import id


def track_ingest(file_id):

    engine = get_engine()
    metadata = sql.MetaData()

    from nbproject import meta

    from lamindb._configuration import user_id

    interface_id = meta.id

    track_do = sql.Table(
        "track_do",
        metadata,
        sql.Column("id", sql.String, primary_key=True, default=id.id_track),
        sql.Column("time", sql.DateTime, default=sql.sql.func.now()),
        sql.Column("user", sql.String, default=user_id),
        sql.Column("interface", sql.String, default=interface_id),
        autoload_with=engine,
    )

    with engine.begin() as conn:
        stmt = sql.insert(track_do).values(type="ingest", file=file_id)
        result = conn.execute(stmt)

    return result.inserted_primary_key[0]


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
