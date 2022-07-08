import shutil
from pathlib import Path
from typing import Union

from cloudpathlib import CloudPath

from lamindb import _setup


def store_file(filepath: Union[str, Path], filekey: str):
    """Store arbitrary file."""
    storage_path = _setup._settings.storage_filepath(filekey)
    if isinstance(storage_path, CloudPath):
        storage_path.upload_from(filepath)
    else:
        shutil.copyfile(filepath, storage_path)
