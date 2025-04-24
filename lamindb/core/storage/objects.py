from __future__ import annotations

from pathlib import PurePosixPath
from typing import TYPE_CHECKING, TypeAlias

from anndata import AnnData
from pandas import DataFrame

from lamindb.core._compat import (
    with_package_obj,
)
from lamindb.core.types import ScverseDataStructures

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr

SupportedDataTypes: TypeAlias = DataFrame | ScverseDataStructures


def infer_suffix(dmem: SupportedDataTypes, format: str | None = None):
    """Infer LaminDB storage file suffix from a data object."""
    if isinstance(dmem, AnnData):
        if format is not None:
            # should be `.h5ad`, `.`zarr`, or `.anndata.zarr`
            if format not in {"h5ad", "zarr", "anndata.zarr"}:
                raise ValueError(
                    "Error when specifying AnnData storage format, it should be"
                    f" 'h5ad', 'zarr', not '{format}'. Check 'format'"
                    " or the suffix of 'key'."
                )
            return "." + format
        return ".h5ad"

    if isinstance(dmem, DataFrame):
        if format == ".csv":
            return ".csv"
        return ".parquet"

    if with_package_obj(
        dmem,
        "MuData",
        "mudata",
        lambda obj: True,  # Just checking type, not calling any method
    )[0]:
        return ".h5mu"

    has_spatialdata, spatialdata_suffix = with_package_obj(
        dmem,
        "SpatialData",
        "spatialdata",
        lambda obj: (
            format
            if format is not None and format in {"spatialdata.zarr", "zarr"}
            else ".zarr"
            if format is None
            else (_ for _ in ()).throw(
                ValueError(
                    "Error when specifying SpatialData storage format, it should be"
                    f" 'zarr', 'spatialdata.zarr', not '{format}'. Check 'format'"
                    " or the suffix of 'key'."
                )
            )
        ),
    )
    if has_spatialdata:
        return spatialdata_suffix
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
        if filepath.suffix == ".csv":
            dmem.to_csv(filepath)
            return
        dmem.to_parquet(filepath)
        return

    if with_package_obj(dmem, "MuData", "mudata", lambda obj: obj.write(filepath))[0]:
        return

    if with_package_obj(
        dmem,
        "SpatialData",
        "spatialdata",
        lambda obj: obj.write(filepath, overwrite=True),
    )[0]:
        return

    raise NotImplementedError
