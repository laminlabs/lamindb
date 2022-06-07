import shutil
from pathlib import Path

import nbproject

from . import db


def ingest(filepath):
    from lamindb._configuration import storage_root

    storage_root = Path(storage_root)

    filepath = Path(filepath)
    filename = filepath.name

    file_id = db.insert.file(filename, nbproject.meta.uid)

    storage_name = str(filepath.stem) + f"--lndb-{file_id}" + str(filepath.suffix)

    storage_path = storage_root / storage_name

    shutil.copyfile(filepath, storage_path)

    print(f"ingested file {file_id} from notebook {nbproject.meta.uid}")
