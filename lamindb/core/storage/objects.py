from __future__ import annotations

from pathlib import PurePosixPath
from typing import TYPE_CHECKING, TypeAlias, TypeVar

from anndata import AnnData
from pandas import DataFrame

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr

SpatialData = TypeVar("SpatialData")
MuData = TypeVar("MuData")

SupportedDataTypes: TypeAlias = AnnData | DataFrame | MuData | SpatialData


def is_package_installed(package_name):
    import importlib.util

    spec = importlib.util.find_spec(package_name)
    return spec is not None


def infer_suffix(dmem: SupportedDataTypes, format: str | None = None):
    """Infer LaminDB storage file suffix from a data object."""
    if isinstance(dmem, AnnData):
        if format is not None:
            if format not in {"h5ad", "zarr", "anndata.zarr"}:
                raise ValueError(
                    "Error when specifying AnnData storage format, it should be"
                    f" 'h5ad', 'zarr', not '{format}'. Check 'format'"
                    " or the suffix of 'key'."
                )
            return "." + format
        return ".h5ad"

    if isinstance(dmem, DataFrame):
        return ".parquet"

    if is_package_installed("mudata"):
        from mudata import MuData

        if isinstance(dmem, MuData):
            return ".h5mu"

    if is_package_installed("spatialdata"):
        from spatialdata import SpatialData

        if isinstance(dmem, SpatialData):
            if format is not None:
                if format not in {"spatialdata.zarr"}:
                    raise ValueError(
                        "Error when specifying SpatialData storage format, it should be"
                        f" 'zarr', 'spatialdata.zarr', not '{format}'. Check 'format'"
                        " or the suffix of 'key'."
                    )
                return "." + format
            return ".zarr"
    else:
        raise NotImplementedError


def write_to_disk(dmem: SupportedDataTypes, filepath: UPathStr) -> None:
    """Writes the passed in memory data to disk to a specified path."""
    if isinstance(dmem, AnnData):
        suffix = PurePosixPath(filepath).suffix
        if suffix == ".h5ad":
            dmem.write_h5ad(filepath)
            return
        elif suffix == ".zarr":
            dmem.write_zarr(filepath)
            return
        else:
            raise NotImplementedError

    if isinstance(dmem, DataFrame):
        dmem.to_parquet(filepath)
        return

    if is_package_installed("mudata"):
        from mudata import MuData

        if isinstance(dmem, MuData):
            dmem.write(filepath)
            return

    if is_package_installed("spatialdata"):
        from spatialdata import SpatialData

        if isinstance(dmem, SpatialData):
            dmem.write(filepath, overwrite=True)
            return
    else:
        raise NotImplementedError
