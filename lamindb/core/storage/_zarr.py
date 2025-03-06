from __future__ import annotations

import warnings
from typing import TYPE_CHECKING, Literal

import scipy.sparse as sparse
import zarr
from anndata import __version__ as anndata_version
from anndata._io.specs import write_elem
from fsspec.implementations.local import LocalFileSystem
from lamin_utils import logger
from lamindb_setup.core.upath import S3FSMap, create_mapper, infer_filesystem
from packaging import version

from lamindb.core._compat import with_package

from ._anndata_sizes import _size_elem, _size_raw, size_adata

if version.parse(anndata_version) < version.parse("0.11.0"):
    from anndata._io import read_zarr as read_anndata_zarr
else:
    from anndata.io import read_zarr as read_anndata_zarr


if TYPE_CHECKING:
    from anndata import AnnData
    from fsspec import FSMap
    from lamindb_setup.core.types import UPathStr

    from lamindb.core.types import ScverseDataStructures


def create_zarr_open_obj(
    storepath: UPathStr, *, check: bool = True
) -> str | S3FSMap | FSMap:
    """Creates the correct object that can be used to open a zarr file depending on local or remote location."""
    fs, storepath_str = infer_filesystem(storepath)

    if isinstance(fs, LocalFileSystem):
        open_obj = storepath_str
    else:
        open_obj = create_mapper(fs, storepath_str, check=check)

    return open_obj


def _identify_zarr_type_from_storage(
    storage: zarr.Group,
) -> Literal["anndata", "mudata", "spatialdata", "unknown"]:
    """Internal helper to identify zarr type from an open storage object."""
    try:
        if storage.attrs.get("encoding-type", "") == "anndata":
            return "anndata"
        elif storage.attrs.get("encoding-type", "") == "MuData":
            return "mudata"
        elif "spatialdata_attrs" in storage.attrs:
            return "spatialdata"
    except Exception as error:
        logger.warning(f"an exception occurred {error}")
    return "unknown"


def identify_zarr_type(
    storepath: UPathStr, *, check: bool = True
) -> Literal["anndata", "mudata", "spatialdata", "unknown"]:
    """Identify whether a zarr store is AnnData, SpatialData, or unknown type."""
    # we can add these cheap suffix-based-checks later
    # also need to check whether the .spatialdata.zarr suffix
    # actually becomes a "standard"; currently we don't recognize it
    # unlike ".anndata.zarr" in VALID_SUFFIXES
    # suffixes = UPath(storepath).suffixes
    # if ".spatialdata" in suffixes:
    #     return "spatialdata"
    # elif ".anndata" in suffixes:
    #     return "anndata"

    open_obj = create_zarr_open_obj(storepath, check=check)
    try:
        storage = zarr.open(open_obj, mode="r")
        return _identify_zarr_type_from_storage(storage)
    except Exception as error:
        logger.warning(
            f"an exception occured while trying to open the zarr store\n {error}"
        )
    return "unknown"


def load_zarr(
    storepath: UPathStr,
    expected_type: Literal["anndata", "mudata", "spatialdata"] = None,
) -> ScverseDataStructures:
    """Loads a zarr store and returns the corresponding scverse data structure.

    Args:
        storepath: Path to the zarr store
        expected_type: If provided, ensures the zarr store is of this type ("anndata", "mudata", "spatialdata")
                       and raises ValueError if it's not
    """
    open_obj = create_zarr_open_obj(storepath, check=True)

    # Open the storage once
    try:
        storage = zarr.open(open_obj, mode="r")
    except Exception as error:
        raise ValueError(f"Could not open zarr store: {error}") from None

    actual_type = _identify_zarr_type_from_storage(storage)
    if expected_type is not None and actual_type != expected_type:
        raise ValueError(
            f"Expected zarr store of type '{expected_type}', but found '{actual_type}'"
        )

    match actual_type:
        case "anndata":
            scverse_obj = read_anndata_zarr(open_obj)
        case "mudata":
            scverse_obj = with_package("mudata", lambda mod: mod.read_zarr(open_obj))
        case "spatialdata":
            scverse_obj = with_package(
                "spatialdata", lambda mod: mod.read_zarr(open_obj)
            )
        case "unknown" | _:
            raise ValueError(
                "Unable to determine zarr store format and therefore cannot load Artifact."
            )
    return scverse_obj


def write_adata_zarr(
    adata: AnnData, storepath: UPathStr, callback=None, chunks=None, **dataset_kwargs
) -> None:
    fs, storepath_str = infer_filesystem(storepath)
    store = create_mapper(fs, storepath_str, create=True)

    f = zarr.open(store, mode="w")

    adata.strings_to_categoricals()
    if adata.raw is not None:
        adata.strings_to_categoricals(adata.raw.var)

    f.attrs.setdefault("encoding-type", "anndata")
    f.attrs.setdefault("encoding-version", "0.1.0")

    adata_size = None
    cumulative_val = 0

    def _report_progress(key_write: str | None = None):
        nonlocal adata_size
        nonlocal cumulative_val

        if callback is None:
            return None
        if adata_size is None:
            adata_size = size_adata(adata)
        if key_write is None:
            # begin or finish
            if cumulative_val < adata_size:
                callback(adata_size, adata_size if cumulative_val > 0 else 0)
            return None

        elem = getattr(adata, key_write, None)
        if elem is None:
            return None
        elem_size = _size_raw(elem) if key_write == "raw" else _size_elem(elem)
        if elem_size == 0:
            return None

        cumulative_val += elem_size
        callback(adata_size, cumulative_val)

    def _write_elem_cb(f, k, elem, dataset_kwargs):
        write_elem(f, k, elem, dataset_kwargs=dataset_kwargs)
        _report_progress(k)

    _report_progress(None)
    with warnings.catch_warnings():
        warnings.filterwarnings("ignore", category=UserWarning, module="zarr")

        if chunks is not None and not isinstance(adata.X, sparse.spmatrix):
            _write_elem_cb(
                f,
                "X",
                adata.X,
                dataset_kwargs=dict(chunks=chunks, **dataset_kwargs),
            )
        else:
            _write_elem_cb(f, "X", adata.X, dataset_kwargs=dataset_kwargs)
        for elem in ("obs", "var"):
            _write_elem_cb(f, elem, getattr(adata, elem), dataset_kwargs=dataset_kwargs)
        for elem in ("obsm", "varm", "obsp", "varp", "layers", "uns"):
            _write_elem_cb(
                f, elem, dict(getattr(adata, elem)), dataset_kwargs=dataset_kwargs
            )
        _write_elem_cb(f, "raw", adata.raw, dataset_kwargs=dataset_kwargs)
    # todo: fix size less than total at the end
    _report_progress(None)
