from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Literal

from anndata._io.specs.registry import get_spec
from lnschema_core import Artifact

from ._anndata_accessor import AnnDataAccessor, StorageType, registry
from .paths import filepath_from_artifact

if TYPE_CHECKING:
    from fsspec.core import OpenFile
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from upath import UPath


def _open_tiledbsoma(
    filepath: UPath, mode: Literal["r", "w"] = "r"
) -> SOMACollection | SOMAExperiment:
    try:
        import tiledbsoma as soma
    except ImportError as e:
        raise ImportError("Please install tiledbsoma: pip install tiledbsoma") from e
    filepath_str = filepath.as_posix()
    if filepath.protocol == "s3":
        from lamindb_setup.core._settings_storage import get_storage_region

        region = get_storage_region(filepath_str)
        tiledb_config = {"vfs.s3.region": region}
        storage_options = filepath.storage_options
        if "key" in storage_options:
            tiledb_config["vfs.s3.aws_access_key_id"] = storage_options["key"]
        if "secret" in storage_options:
            tiledb_config["vfs.s3.aws_secret_access_key"] = storage_options["secret"]
        if "token" in storage_options:
            tiledb_config["vfs.s3.aws_session_token"] = storage_options["token"]
        ctx = soma.SOMATileDBContext(tiledb_config=tiledb_config)
        # this is a strange bug
        # for some reason iterdir futher gives incorrect results
        # if cache is not invalidated
        # instead of obs and ms it gives ms and ms in the list of names
        filepath.fs.invalidate_cache()
    else:
        ctx = None

    soma_objects = [obj.name for obj in filepath.iterdir()]
    if "obs" in soma_objects and "ms" in soma_objects:
        SOMAType = soma.Experiment
    else:
        SOMAType = soma.Collection
    return SOMAType.open(filepath_str, mode=mode, context=ctx)


@dataclass
class BackedAccessor:
    """h5py.File or zarr.Group accessor."""

    connection: OpenFile
    """The connection."""
    storage: StorageType
    """The storage access."""


def backed_access(
    artifact_or_filepath: Artifact | UPath,
    mode: str = "r",
    using_key: str | None = None,
) -> AnnDataAccessor | BackedAccessor | SOMACollection | SOMAExperiment:
    if isinstance(artifact_or_filepath, Artifact):
        filepath = filepath_from_artifact(artifact_or_filepath, using_key=using_key)
    else:
        filepath = artifact_or_filepath
    name = filepath.name
    suffix = filepath.suffix

    if name == "soma" or suffix == ".tiledbsoma":
        if mode not in {"r", "w"}:
            raise ValueError("`mode` should be either 'r' or 'w' for tiledbsoma.")
        return _open_tiledbsoma(filepath, mode=mode)  # type: ignore
    elif suffix in {".h5", ".hdf5", ".h5ad"}:
        conn, storage = registry.open("h5py", filepath, mode=mode)
    elif suffix == ".zarr":
        conn, storage = registry.open("zarr", filepath, mode=mode)
    else:
        raise ValueError(
            "object should have .h5, .hdf5, .h5ad, .zarr, .tiledbsoma suffix, not"
            f" {suffix}."
        )

    is_anndata = suffix == ".h5ad" or get_spec(storage).encoding_type == "anndata"
    if is_anndata:
        if mode != "r":
            raise ValueError("Can only access `AnnData` with mode='r'.")
        return AnnDataAccessor(conn, storage, name)
    else:
        return BackedAccessor(conn, storage)
