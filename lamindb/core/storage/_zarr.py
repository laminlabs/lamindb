from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import zarr
from anndata import __version__ as anndata_version
from lamin_utils import logger
from lamindb_setup.core.upath import LocalPathClasses, S3FSMap, UPath, create_mapper
from packaging import version

from lamindb.core._compat import with_package

if version.parse(anndata_version) < version.parse("0.11.0"):
    from anndata._io import read_zarr as read_anndata_zarr
else:
    from anndata.io import read_zarr as read_anndata_zarr

if version.parse(zarr.__version__) >= version.parse("3.0.0a0"):
    IS_ZARR_V3 = True
    from zarr.abc.store import Store
else:
    IS_ZARR_V3 = False
    from zarr.storage import Store  # noqa

if TYPE_CHECKING:
    from fsspec import FSMap
    from lamindb_setup.core.types import UPathStr

    from lamindb.core.types import ScverseDataStructures


def get_zarr_store(
    path: UPathStr, *, check: bool = False, create: bool = False
) -> str | S3FSMap | FSMap | Store:
    """Creates the correct object that can be used to open a zarr file depending on local or remote location."""
    storepath, storepath_str = UPath(path), str(path)
    if isinstance(storepath, LocalPathClasses):
        store = storepath_str
    elif IS_ZARR_V3:
        store = zarr.storage.FsspecStore.from_upath(UPath(storepath, asynchronous=True))
    else:
        store = create_mapper(storepath.fs, storepath_str, check=check, create=create)

    return store


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
    suffixes = UPath(storepath).suffixes
    if ".anndata" in suffixes:
        return "anndata"
    elif ".mudata" in suffixes:
        return "mudata"
    elif ".spatialdata" in suffixes:
        return "spatialdata"

    store = get_zarr_store(storepath, check=check)
    try:
        storage = zarr.open(store, mode="r")
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
    store = get_zarr_store(storepath, check=True)
    # Open the storage once
    try:
        storage = zarr.open(store, mode="r")
    except Exception as error:
        raise ValueError(f"Could not open zarr store: {error}") from None

    actual_type = _identify_zarr_type_from_storage(storage)
    if expected_type is not None and actual_type != expected_type:
        raise ValueError(
            f"Expected zarr store of type '{expected_type}', but found '{actual_type}'"
        )

    match actual_type:
        case "anndata":
            scverse_obj = read_anndata_zarr(store)
        case "mudata":
            scverse_obj = with_package("mudata", lambda mod: mod.read_zarr(store))
        case "spatialdata":
            scverse_obj = with_package("spatialdata", lambda mod: mod.read_zarr(store))
        case "unknown" | _:
            raise ValueError(
                "Unable to determine zarr store format and therefore cannot load Artifact."
            )
    return scverse_obj
