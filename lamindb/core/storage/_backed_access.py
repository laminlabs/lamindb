from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TYPE_CHECKING, Any, Callable, Literal

from anndata._io.specs.registry import get_spec

from ._anndata_accessor import AnnDataAccessor, StorageType, registry
from ._polars_lazy_df import POLARS_SUFFIXES, _open_polars_lazy_df
from ._pyarrow_dataset import PYARROW_SUFFIXES, _open_pyarrow_dataset
from ._tiledbsoma import _open_tiledbsoma
from .paths import filepath_from_artifact

if TYPE_CHECKING:
    from collections.abc import Iterator

    from fsspec.core import OpenFile
    from polars import LazyFrame as PolarsLazyFrame
    from pyarrow.dataset import Dataset as PyArrowDataset
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
    from tiledbsoma import Measurement as SOMAMeasurement
    from upath import UPath

    from lamindb.models.artifact import Artifact


# this dynamically creates a subclass of a context manager class
# and reassigns it to an instance of the superclass
# so that the instance calls finalize on close or exit
def _track_writes_factory(obj: Any, finalize: Callable):
    closed: bool = False

    tracked_class = obj.__class__
    type_dict = {"__doc__": tracked_class.__doc__}
    if hasattr(tracked_class, "__slots__"):
        type_dict["__slots__"] = ()
    if hasattr(tracked_class, "__exit__"):

        def __exit__(self, exc_type, exc_val, exc_tb):
            nonlocal closed
            tracked_class.__exit__(self, exc_type, exc_val, exc_tb)
            if not closed:
                finalize()
                closed = True

        type_dict["__exit__"] = __exit__
    if hasattr(tracked_class, "close"):

        def close(self, *args, **kwargs):
            nonlocal closed
            tracked_class.close(self, *args, **kwargs)
            if not closed:
                finalize()
                closed = True

        type_dict["close"] = close

    Track = type(tracked_class.__name__ + "Track", (tracked_class,), type_dict)
    obj.__class__ = Track
    return obj


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
    engine: Literal["pyarrow", "polars"] = "pyarrow",
    using_key: str | None = None,
    **kwargs,
) -> (
    AnnDataAccessor
    | BackedAccessor
    | SOMACollection
    | SOMAExperiment
    | SOMAMeasurement
    | PyArrowDataset
    | Iterator[PolarsLazyFrame]
):
    from lamindb.models import Artifact

    if isinstance(artifact_or_filepath, Artifact):
        objectpath, _ = filepath_from_artifact(
            artifact_or_filepath, using_key=using_key
        )
    else:
        objectpath = artifact_or_filepath
    name = objectpath.name
    # ignore .gz, only check the real suffix
    suffixes = objectpath.suffixes
    suffix = (
        suffixes[-2] if len(suffixes) > 1 and ".gz" in suffixes else objectpath.suffix
    )

    if name == "soma" or suffix == ".tiledbsoma":
        if mode not in {"r", "w"}:
            raise ValueError("`mode` should be either 'r' or 'w' for tiledbsoma.")
        return _open_tiledbsoma(objectpath, mode=mode, **kwargs)  # type: ignore
    elif suffix in {".h5", ".hdf5", ".h5ad"}:
        conn, storage = registry.open("h5py", objectpath, mode=mode, **kwargs)
    elif suffix == ".zarr":
        conn, storage = registry.open("zarr", objectpath, mode=mode, **kwargs)
    elif len(df_suffixes := _flat_suffixes(objectpath)) == 1 and (
        df_suffix := df_suffixes.pop()
    ) in set(PYARROW_SUFFIXES).union(POLARS_SUFFIXES):
        return _open_dataframe(objectpath, df_suffix, engine, **kwargs)
    else:
        raise ValueError(
            "The object should have .h5, .hdf5, .h5ad, .zarr, .tiledbsoma suffix "
            f"be compatible with pyarrow.dataset.dataset or polars.scan_* functions, "
            f"instead of being {suffix} object."
        )

    is_anndata = suffix == ".h5ad" or get_spec(storage).encoding_type == "anndata"
    if is_anndata:
        if mode != "r":
            raise ValueError("Can only access `AnnData` with mode='r'.")
        return AnnDataAccessor(conn, storage, name)
    else:
        return BackedAccessor(conn, storage)


def _flat_suffixes(paths: UPath | list[UPath]) -> set[str]:
    # it is assumed here that the paths exist
    # we don't check here that the filesystem is the same
    # but this is a requirement for pyarrow.dataset.dataset
    path_list = []
    if isinstance(paths, Path):
        paths = [paths]
    for path in paths:
        # assume http is always a file
        if getattr(path, "protocol", None) not in {"http", "https"} and path.is_dir():
            path_list += [p for p in path.rglob("*") if p.suffix != ""]
        else:
            path_list.append(path)

    suffixes = set()
    for path in path_list:
        path_suffixes = path.suffixes
        # this doesn't work for externally gzipped files, REMOVE LATER
        path_suffix = (
            path_suffixes[-2]
            if len(path_suffixes) > 1 and ".gz" in path_suffixes
            else path.suffix
        )
        suffixes.add(path_suffix)
    return suffixes


def _open_dataframe(
    paths: UPath | list[UPath],
    suffix: str | None = None,
    engine: Literal["pyarrow", "polars"] = "pyarrow",
    **kwargs,
) -> PyArrowDataset | Iterator[PolarsLazyFrame]:
    df_suffix: str
    if suffix is None:
        df_suffixes = _flat_suffixes(paths)
        if len(df_suffixes) > 1:
            raise ValueError(
                f"The artifacts in the collection have different file formats: {', '.join(df_suffixes)}.\n"
                "It is not possible to open such stores with pyarrow or polars."
            )
        df_suffix = df_suffixes.pop()
    else:
        df_suffix = suffix

    if engine == "pyarrow":
        if df_suffix not in PYARROW_SUFFIXES:
            raise ValueError(
                f"{df_suffix} files are not supported by pyarrow, "
                f"they should have one of these formats: {', '.join(PYARROW_SUFFIXES)}."
            )
        # this checks that the filesystem is the same for all paths
        # this is a requirement of pyarrow.dataset.dataset
        if not isinstance(paths, Path):  # is a list then
            fs = getattr(paths[0], "fs", None)
            for path in paths[1:]:
                # this assumes that the filesystems are cached by fsspec
                if getattr(path, "fs", None) is not fs:
                    raise ValueError(
                        "The collection has artifacts with different filesystems, "
                        "this is not supported by pyarrow."
                    )
        dataframe = _open_pyarrow_dataset(paths, **kwargs)
    elif engine == "polars":
        if df_suffix not in POLARS_SUFFIXES:
            raise ValueError(
                f"{df_suffix} files are not supported by polars, "
                f"they should have one of these formats: {', '.join(POLARS_SUFFIXES)}."
            )
        dataframe = _open_polars_lazy_df(paths, **kwargs)
    else:
        raise ValueError(
            f"Unknown engine: {engine}. It should be 'pyarrow' or 'polars'."
        )

    return dataframe
