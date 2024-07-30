from __future__ import annotations

import builtins
import re
import shutil
from pathlib import Path
from typing import TYPE_CHECKING

import anndata as ad
import pandas as pd
from lamin_utils import logger
from lamindb_setup import settings as setup_settings
from lamindb_setup.core import StorageSettings
from lamindb_setup.core.upath import (
    LocalPathClasses,
    UPath,
    create_path,
    infer_filesystem,
)
from lnschema_core.models import Artifact, Storage

from lamindb.core._settings import settings

if TYPE_CHECKING:
    import mudata as md
    from lamindb_setup.core.types import UPathStr

try:
    from ._zarr import read_adata_zarr
except ImportError:

    def read_adata_zarr(storepath):  # type: ignore
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
    assert isinstance(suffix, str)  # noqa: S101 Suffix cannot be None.
    if is_dir:
        uid_storage = uid[:16]  # 16 chars, leave 4 chars for versioning
    else:
        uid_storage = uid
    storage_key = f"{AUTO_KEY_PREFIX}{uid_storage}{suffix}"
    return storage_key


def check_path_is_child_of_root(path: Path | UPath, root: Path | UPath | None) -> bool:
    # str is needed to eliminate UPath storage_options
    # from the equality checks below
    path = UPath(str(path))
    root = UPath(str(root))
    return root.resolve() in path.resolve().parents


def attempt_accessing_path(
    artifact: Artifact,
    storage_key: str,
    using_key: str | None = None,
    access_token: str | None = None,
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
                settings.storage.root, access_token=access_token
            )
    else:
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
def filepath_from_artifact(artifact: Artifact, using_key: str | None = None):
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


def store_file_or_folder(
    local_path: UPathStr, storage_path: UPath, print_progress: bool = True
) -> None:
    """Store file or folder (localpath) at storagepath."""
    local_path = Path(local_path)
    if not isinstance(storage_path, LocalPathClasses):
        # this uploads files and directories
        create_folder = False if local_path.is_dir() else None
        storage_path.upload_from(
            local_path, create_folder=create_folder, print_progress=print_progress
        )
    else:  # storage path is local
        storage_path.parent.mkdir(parents=True, exist_ok=True)
        if local_path.is_file():
            try:
                shutil.copyfile(local_path, storage_path)
            except shutil.SameFileError:
                pass
        else:
            if storage_path.exists():
                shutil.rmtree(storage_path)
            shutil.copytree(local_path, storage_path)


def delete_storage_using_key(
    artifact: Artifact, storage_key: str, using_key: str | None
):
    filepath = attempt_accessing_path(artifact, storage_key, using_key=using_key)
    delete_storage(filepath)


def delete_storage(
    storagepath: Path, raise_file_not_found_error: bool = True
) -> None | str:
    """Delete arbitrary artifact."""
    if storagepath.is_file():
        storagepath.unlink()
    elif storagepath.is_dir():
        if isinstance(storagepath, LocalPathClasses) or not isinstance(
            storagepath, UPath
        ):
            shutil.rmtree(storagepath)
        else:
            storagepath.rmdir()
    elif raise_file_not_found_error:
        raise FileNotFoundError(f"{storagepath} is not an existing path!")
    else:
        logger.warning(f"{storagepath} is not an existing path!")
    return None


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


def read_mdata_h5mu(filepath: UPathStr, **kwargs) -> md.MuData:
    import mudata as md

    path_sanitized = Path(filepath)
    return md.read_h5mu(path_sanitized, **kwargs)


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

    if filepath.suffix not in {".h5ad", ".zarr"}:
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
        ".html": load_html,
        ".json": load_json,
        ".h5mu": read_mdata_h5mu,
    }

    reader = READER_FUNCS.get(filepath.suffix)
    if reader is None:
        return filepath
    else:
        return reader(filepath, **kwargs)
