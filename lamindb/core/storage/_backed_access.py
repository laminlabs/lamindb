from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable

from anndata._io.specs.registry import get_spec

from ._anndata_accessor import AnnDataAccessor, StorageType, registry
from ._pyarrow_dataset import _is_pyarrow_dataset, _open_pyarrow_dataset
from ._tiledbsoma import _open_tiledbsoma
from .paths import filepath_from_artifact

if TYPE_CHECKING:
    from fsspec.core import OpenFile
    from pyarrow.dataset import Dataset as PyArrowDataset
    from tiledbsoma import Collection as SOMACollection
    from tiledbsoma import Experiment as SOMAExperiment
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
    using_key: str | None = None,
    **kwargs,
) -> (
    AnnDataAccessor | BackedAccessor | SOMACollection | SOMAExperiment | PyArrowDataset
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
    elif _is_pyarrow_dataset(objectpath):
        return _open_pyarrow_dataset(objectpath, **kwargs)
    else:
        raise ValueError(
            "The object should have .h5, .hdf5, .h5ad, .zarr, .tiledbsoma suffix "
            f"or be compatible with pyarrow.dataset.dataset, instead of being {suffix} object."
        )

    is_anndata = suffix == ".h5ad" or get_spec(storage).encoding_type == "anndata"
    if is_anndata:
        if mode != "r":
            raise ValueError("Can only access `AnnData` with mode='r'.")
        return AnnDataAccessor(conn, storage, name)
    else:
        return BackedAccessor(conn, storage)
