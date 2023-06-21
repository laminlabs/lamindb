import os
import shutil
from pathlib import Path
from typing import Union

import anndata as ad
import fsspec
import pandas as pd
from lamindb_setup import settings
from lamindb_setup.dev.upath import UPath, infer_filesystem
from lnschema_core.models import File

from lamindb._file_access import attempt_accessing_path

from ._zarr import read_adata_zarr


def read_adata_h5ad(filepath, **kwargs) -> ad.AnnData:
    fs, filepath = infer_filesystem(filepath)

    with fs.open(filepath, mode="rb") as file:
        adata = ad.read_h5ad(file, backed=False, **kwargs)
        return adata


def print_hook(size, value, **kwargs):
    progress = value / size
    out = f"... uploading {Path(kwargs['filepath']).name}: {min(progress, 1.):4.2f}"
    if progress >= 1:
        out += "\n"
    if "NBPRJ_TEST_NBPATH" not in os.environ:
        print(out, end="\r")


class ProgressCallback(fsspec.callbacks.Callback):
    def branch(self, path_1, path_2, kwargs):
        kwargs["callback"] = fsspec.callbacks.Callback(
            hooks=dict(print_hook=print_hook), filepath=path_1
        )

    def call(self, *args, **kwargs):
        return None


def store_object(localpath: Union[str, Path], storagekey: str) -> float:
    """Store arbitrary file to configured storage location.

    Returns size in bytes.
    """
    storagepath = settings.instance.storage.key_to_filepath(storagekey)
    localpath = Path(localpath)

    if localpath.is_file():
        size = localpath.stat().st_size
    else:
        size = sum(f.stat().st_size for f in localpath.rglob("*") if f.is_file())

    if isinstance(storagepath, UPath):
        if localpath.suffix != ".zarr":
            cb = ProgressCallback()
        else:
            # todo: make proper progress bar for zarr
            cb = fsspec.callbacks.NoOpCallback()
        storagepath.upload_from(localpath, recursive=True, callback=cb)
    else:  # storage path is local
        storagepath.parent.mkdir(parents=True, exist_ok=True)
        if localpath.is_file():
            try:
                shutil.copyfile(localpath, storagepath)
            except shutil.SameFileError:
                pass
        else:
            if storagepath.exists():
                shutil.rmtree(storagepath)
            shutil.copytree(localpath, storagepath)
    return float(size)  # because this is how we store in the db


def delete_storage_using_key(file: File, storage_key: str):
    filepath = attempt_accessing_path(file, storage_key)
    delete_storage(filepath)


def delete_storage(storagepath: Union[Path, UPath]):
    """Delete arbitrary file."""
    if storagepath.is_file():
        storagepath.unlink()
    elif storagepath.is_dir():
        if isinstance(storagepath, UPath):
            storagepath.rmdir()
        else:
            shutil.rmtree(storagepath)
    else:
        raise FileNotFoundError(f"{storagepath} is not an existing path!")


def read_fcs(*args, **kwargs):
    try:
        import readfcs
    except ImportError:
        raise ImportError("Please install readfcs: pip install readfcs")
    return readfcs.read(*args, **kwargs)


def load_to_memory(filepath: Union[str, Path, UPath], stream: bool = False):
    """Load a file into memory.

    Returns the filepath if no in-memory form is found.
    """
    if isinstance(filepath, str):
        filepath = Path(filepath)

    if filepath.suffix in (".zarr", ".zrad"):
        stream = True
    elif filepath.suffix != ".h5ad":
        stream = False

    if not stream:
        # caching happens here if filename is a UPath
        # todo: make it safe when filepath is just Path
        filepath = settings.instance.storage.cloud_to_local(filepath)

    READER_FUNCS = {
        ".csv": pd.read_csv,
        ".h5ad": read_adata_h5ad,
        ".parquet": pd.read_parquet,
        ".fcs": read_fcs,
        ".zarr": read_adata_zarr,
        ".zrad": read_adata_zarr,
    }

    reader = READER_FUNCS.get(filepath.suffix)
    if reader is None:
        return filepath
    else:
        return reader(filepath)
