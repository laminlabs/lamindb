import shutil
from pathlib import Path
from typing import Union

import pandas as pd
import readfcs
from cloudpathlib import CloudPath
from lndb_setup import settings

from ._h5ad import read_adata_h5ad
from ._zarr import read_adata_zarr

READER_FUNCS = {
    ".csv": pd.read_csv,
    ".h5ad": read_adata_h5ad,
    ".feather": pd.read_feather,
    ".fcs": readfcs.read,
    ".zarr": read_adata_zarr,
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


def load_to_memory(filepath: Union[str, Path], stream: bool = False):
    """Load a file into memory.

    Returns the filepath if no in-memory form is found.
    """
    if isinstance(filepath, str):
        filepath = Path(filepath)

    if filepath.suffix == ".zarr":
        stream = True
    elif filepath.suffix != ".h5ad":
        stream = False

    if not stream:
        # caching happens here if filename is a CloudPath
        filepath = Path(filepath)

    reader = READER_FUNCS.get(filepath.suffix)
    if reader is None:
        return filepath
    else:
        return reader(filepath)
