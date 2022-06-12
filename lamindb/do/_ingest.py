import shutil
from pathlib import Path


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

    from ._admin.db import insert

    file_id = insert.file(filename)

    storage_name = str(filepath.stem) + f"--lndb-{file_id}" + str(filepath.suffix)

    storage_path = storage_root / storage_name

    shutil.copyfile(filepath, storage_path)

    import nbproject

    from lamindb._configuration import user_id

    print(f"added file {file_id} from source {nbproject.meta.id} by user {user_id}")
