import shutil
from pathlib import Path
from typing import Union

import fsspec
import pandas as pd
import readfcs
from cloudpathlib import CloudPath
from lndb_setup import settings

from ._h5ad import read_adata_h5ad
from ._zarr import read_adata_zarr

READER_FUNCS = {
    ".csv": pd.read_csv,
    ".h5ad": read_adata_h5ad,
    ".parquet": pd.read_parquet,
    ".fcs": readfcs.read,
    ".zarr": read_adata_zarr,
}

fsspec_filesystem = None


def print_hook(size, value, **kwargs):
    progress = value / size
    out = f"Writing {kwargs['filepath']}: {min(progress, 1.):4.2f}"
    if progress >= 1:
        out += "\n"
    print(out, end="\r")


class ProgressCallback(fsspec.callbacks.Callback):
    def branch(self, path_1, path_2, kwargs):
        kwargs["callback"] = fsspec.callbacks.Callback(
            hooks=dict(print_hook=print_hook), filepath=path_1
        )

    def call(self, *args, **kwargs):
        return None


def store_file(
    localfile: Union[str, Path], storagekey: str, use_fsspec: bool = False
) -> float:
    """Store arbitrary file.

    Returns size in bytes.
    """
    storagepath = settings.instance.storage.key_to_filepath(storagekey)
    global fsspec_filesystem
    if isinstance(storagepath, CloudPath):
        if use_fsspec:
            if fsspec_filesystem is None:
                fsspec_filesystem = fsspec.filesystem(
                    storagepath.cloud_prefix.replace("://", "")
                )
            fsspec_filesystem.put(
                str(localfile),
                str(storagepath),
                recursive=True,
                callback=ProgressCallback(),
            )
        else:
            storagepath.upload_from(localfile)
    else:
        try:
            shutil.copyfile(localfile, storagepath)
        except shutil.SameFileError:
            pass
    size = Path(localfile).stat().st_size
    return float(size)  # because this is how we store in the db


def delete_storage(storagekey: str):
    """Delete arbitrary file."""
    storagepath = settings.instance.storage.key_to_filepath(storagekey)
    if storagepath.is_file():
        storagepath.unlink()
    else:
        storagepath.rmtree()


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
