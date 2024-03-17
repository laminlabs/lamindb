import builtins
import re
import shutil
from pathlib import Path
from typing import Optional, Union

import anndata as ad
import pandas as pd
from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core import StorageSettings
from lamindb_setup.core.types import UPathStr
from lamindb_setup.core.upath import (
    LocalPathClasses,
    UPath,
    create_path,
    infer_filesystem,
)
from lnschema_core.models import Artifact, Storage

from lamindb.core._settings import settings

try:
    from ._zarr import read_adata_zarr
except ImportError:

    def read_adata_zarr(filepath):  # type: ignore
        raise ImportError("Please install zarr: pip install zarr")


AUTO_KEY_PREFIX = ".lamindb/"
is_run_from_ipython = getattr(builtins, "__IPYTHON__", False)


# add type annotations back asap when re-organizing the module
def auto_storage_key_from_artifact(artifact: Artifact):
    if artifact.key is None or artifact.key_is_virtual:
        is_dir = artifact.n_objects is not None
        return auto_storage_key_from_artifact_uid(artifact.uid, artifact.suffix, is_dir)
    else:
        return artifact.key


def auto_storage_key_from_artifact_uid(uid: str, suffix: str, is_dir: bool) -> str:
    assert isinstance(suffix, str)
    if is_dir:
        uid_storage = uid[:16]  # 16 chars, leave 4 chars for versioning
    else:
        uid_storage = uid
    storage_key = f"{AUTO_KEY_PREFIX}{uid_storage}{suffix}"
    return storage_key


def attempt_accessing_path(
    artifact: Artifact,
    storage_key: str,
    using_key: Optional[str] = None,
    access_token: Optional[str] = None,
):
    # check whether the file is in the default db and whether storage
    # matches default storage
    if (
        artifact._state.db in ("default", None)
        and artifact.storage_id == settings._storage_settings.id
    ):
        if access_token is None:
            storage_settings = settings._storage_settings
        else:
            storage_settings = StorageSettings(
                settings.storage, access_token=access_token
            )
    else:
        logger.debug(
            "artifact.path is slightly slower for files outside default storage"
        )
        if artifact._state.db not in ("default", None) and using_key is None:
            storage = (
                Storage.using(artifact._state.db).filter(id=artifact.storage_id).one()
            )
        else:
            storage = (
                Storage.objects.using(using_key).filter(id=artifact.storage_id).one()
            )
        # find a better way than passing None to instance_settings in the future!
        storage_settings = StorageSettings(storage.root, access_token=access_token)
    path = storage_settings.key_to_filepath(storage_key)
    return path


# add type annotations back asap when re-organizing the module
def filepath_from_artifact(artifact: Artifact, using_key: Optional[str] = None):
    if hasattr(artifact, "_local_filepath") and artifact._local_filepath is not None:
        return artifact._local_filepath.resolve()
    storage_key = auto_storage_key_from_artifact(artifact)
    path = attempt_accessing_path(artifact, storage_key, using_key=using_key)
    return path


def read_adata_h5ad(filepath, **kwargs) -> ad.AnnData:
    fs, filepath = infer_filesystem(filepath)

    with fs.open(filepath, mode="rb") as file:
        adata = ad.read_h5ad(file, backed=False, **kwargs)
        return adata


def store_artifact(localpath: UPathStr, storagepath: UPath) -> None:
    """Store directory or file to configured storage location.

    Returns size in bytes.
    """
    localpath = Path(localpath)
    if not isinstance(storagepath, LocalPathClasses):
        # this uploads files and directories
        storagepath.upload_from(localpath, recursive=True, print_progress=True)
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


def delete_storage_using_key(
    artifact: Artifact, storage_key: str, using_key: Optional[str]
):
    filepath = attempt_accessing_path(artifact, storage_key, using_key=using_key)
    delete_storage(filepath)


def delete_storage(storagepath: Path):
    """Delete arbitrary artifact."""
    if not storagepath.is_relative_to(settings.storage):  # type: ignore
        logger.warning("couldn't delete files outside of default storage")
        return "did-not-delete"
    # only delete files in the default storage
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
        raise ImportError("Please install readfcs: pip install readfcs") from None
    return readfcs.read(*args, **kwargs)


def read_tsv(path: UPathStr, **kwargs) -> pd.DataFrame:
    path_sanitized = Path(path)
    return pd.read_csv(path_sanitized, sep="\t", **kwargs)


def load_html(path: UPathStr):
    if is_run_from_ipython:
        with open(path, encoding="utf-8") as f:
            html_content = f.read()
        # Extract the body content using regular expressions
        body_content = re.findall(
            r"<body(?:.*?)>(?:.*?)</body>", html_content, re.DOTALL
        )
        # Remove any empty body tags
        if body_content:
            body_content = body_content[0]
            body_content = body_content.strip()  # type: ignore
        from IPython.display import HTML, display

        display(HTML(data=body_content))
    else:
        return path


def load_json(path: UPathStr):
    import json

    with open(path) as f:
        data = json.load(f)
    return data


def load_to_memory(filepath: UPathStr, stream: bool = False, **kwargs):
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
        filepath = settings._storage_settings.cloud_to_local(
            filepath, print_progress=True
        )

    READER_FUNCS = {
        ".csv": pd.read_csv,
        ".tsv": read_tsv,
        ".h5ad": read_adata_h5ad,
        ".parquet": pd.read_parquet,
        ".fcs": read_fcs,
        ".zarr": read_adata_zarr,
        ".zrad": read_adata_zarr,
        ".html": load_html,
        ".json": load_json,
    }

    reader = READER_FUNCS.get(filepath.suffix)
    if reader is None:
        return filepath
    else:
        return reader(filepath, **kwargs)
