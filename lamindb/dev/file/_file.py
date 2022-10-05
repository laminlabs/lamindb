import shutil
from pathlib import Path
from typing import Union

import anndata as ad
import pandas as pd
import readfcs
from cloudpathlib import CloudPath
from lndb_setup import settings

READER_FUNCS = {
    ".csv": pd.read_csv,
    ".h5ad": ad.read,
    ".feather": pd.read_feather,
    ".fcs": readfcs.read,
}


def store_file(localfile: Union[str, Path], storagekey: str) -> float:
    """Store arbitrary file.

    Returns size in bytes.
    """
    storagepath = settings.instance.storage.key_to_filepath(storagekey)
    if isinstance(storagepath, CloudPath):
        storagepath.upload_from(localfile)
    else:
        try:
            shutil.copyfile(localfile, storagepath)
        except shutil.SameFileError:
            pass
    size = Path(localfile).stat().st_size
    return float(size)  # because this is how we store in the db


def delete_file(storagekey: str):
    """Delete arbitrary file."""
    storagepath = settings.instance.storage.key_to_filepath(storagekey)
    storagepath.unlink()


def load_to_memory(filepath: Union[str, Path]):
    filepath = Path(filepath)
    reader = READER_FUNCS.get(filepath.suffix)
    if reader is None:
        raise NotImplementedError
    else:
        return reader(filepath)
