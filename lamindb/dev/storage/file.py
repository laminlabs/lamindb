import os
import shutil
from pathlib import Path
from typing import Literal, Union

import anndata as ad
import fsspec
import pandas as pd
from lamin_utils import logger
from lamindb_setup import settings
from lamindb_setup.dev import StorageSettings
from lamindb_setup.dev.upath import (
    LocalPathClasses,
    UPath,
    create_path,
    infer_filesystem,
)
from lnschema_core.models import File, Storage

try:
    from ._zarr import read_adata_zarr
except ImportError:

    def read_adata_zarr(filepath):  # type: ignore
        raise ImportError("Please install zarr: pip install zarr")


AUTO_KEY_PREFIX = ".lamindb/"


# add type annotations back asap when re-organizing the module
def auto_storage_key_from_file(file: File):
    if file.key is None or file.key_is_virtual:
        return auto_storage_key_from_id_suffix(file.uid, file.suffix)
    else:
        return file.key


def auto_storage_key_from_id_suffix(uid: str, suffix: str) -> str:
    assert isinstance(uid, str)
    assert isinstance(suffix, str)
    storage_key = f"{AUTO_KEY_PREFIX}{uid}{suffix}"
    return storage_key


def attempt_accessing_path(file: File, storage_key: str):
    if file.storage_id == settings.storage.id:
        path = settings.storage.key_to_filepath(storage_key)
    else:
        logger.debug("file.path is slightly slower for files outside default storage")
        storage = Storage.filter(id=file.storage_id).one()
        # find a better way than passing None to instance_settings in the future!
        storage_settings = StorageSettings(storage.root)
        path = storage_settings.key_to_filepath(storage_key)
    return path


# add type annotations back asap when re-organizing the module
def filepath_from_file(file: File):
    if hasattr(file, "_local_filepath") and file._local_filepath is not None:
        return file._local_filepath.resolve()
    storage_key = auto_storage_key_from_file(file)
    path = attempt_accessing_path(file, storage_key)
    return path


def read_adata_h5ad(filepath, **kwargs) -> ad.AnnData:
    fs, filepath = infer_filesystem(filepath)

    with fs.open(filepath, mode="rb") as file:
        adata = ad.read_h5ad(file, backed=False, **kwargs)
        return adata


def print_hook(size: int, value: int, **kwargs):
    progress = value / size
    out = (
        f"... {kwargs['action']} {Path(kwargs['filepath']).name}:"
        f" {min(progress, 1.):4.2f}"
    )
    if progress >= 1:
        out += "\n"
    if "NBPRJ_TEST_NBPATH" not in os.environ:
        print(out, end="\r")


class ProgressCallback(fsspec.callbacks.Callback):
    def __init__(self, action: Literal["uploading", "downloading"]):
        super().__init__()
        self.action = action

    def branch(self, path_1, path_2, kwargs):
        kwargs["callback"] = fsspec.callbacks.Callback(
            hooks=dict(print_hook=print_hook), filepath=path_1, action=self.action
        )

    def call(self, *args, **kwargs):
        return None


def store_object(localpath: Union[str, Path, UPath], storagekey: str) -> float:
    """Store arbitrary file to configured storage location.

    Returns size in bytes.
    """
    storagepath: UPath = settings.instance.storage.key_to_filepath(storagekey)
    localpath = Path(localpath)

    if localpath.is_file():
        size = localpath.stat().st_size
    else:
        size = sum(f.stat().st_size for f in localpath.rglob("*") if f.is_file())

    if not isinstance(storagepath, LocalPathClasses):
        if localpath.suffix not in {".zarr", ".zrad"}:
            cb = ProgressCallback("uploading")
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
        if isinstance(storagepath, LocalPathClasses) or not isinstance(
            storagepath, UPath
        ):
            shutil.rmtree(storagepath)
        else:
            storagepath.rmdir()
    else:
        raise FileNotFoundError(f"{storagepath} is not an existing path!")


# tested in lamin-usecases
def read_fcs(*args, **kwargs):
    try:
        import readfcs
    except ImportError:  # pragma: no cover
        raise ImportError("Please install readfcs: pip install readfcs")
    return readfcs.read(*args, **kwargs)


def read_tsv(path: Union[str, Path, UPath]) -> pd.DataFrame:
    path_sanitized = Path(path)
    return pd.read_csv(path_sanitized, sep="\t")


def load_to_memory(filepath: Union[str, Path, UPath], stream: bool = False, **kwargs):
    """Load a file into memory.

    Returns the filepath if no in-memory form is found.
    """
    filepath = create_path(filepath)

    if filepath.suffix in (".zarr", ".zrad"):
        stream = True
    elif filepath.suffix != ".h5ad":
        stream = False

    if not stream:
        # caching happens here if filename is a UPath
        # todo: make it safe when filepath is just Path
        cb = ProgressCallback("downloading")
        filepath = settings.instance.storage.cloud_to_local(filepath, callback=cb)

    READER_FUNCS = {
        ".csv": pd.read_csv,
        ".tsv": read_tsv,
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
        return reader(filepath, **kwargs)
