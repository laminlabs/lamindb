from __future__ import annotations

from pathlib import Path
from typing import TYPE_CHECKING, Any, TypeAlias

from lamindb.core._compat import (
    with_package_obj,
)

if TYPE_CHECKING:
    from lamindb_setup.types import AnyPathStr
    from pandas import DataFrame

    from .types import ScverseDataStructures

    SupportedDataTypes: TypeAlias = DataFrame | ScverseDataStructures
else:
    SupportedDataTypes: TypeAlias = Any


def infer_suffix(
    dmem: SupportedDataTypes, format: str | dict[str, Any] | None = None
) -> str:
    """Infer LaminDB storage file suffix from a data object."""
    has_anndata, anndata_suffix = with_package_obj(
        dmem,
        "AnnData",
        "anndata",
        lambda obj: _infer_anndata_suffix(format),
    )
    if has_anndata:
        return anndata_suffix

    has_dataframe, dataframe_suffix = with_package_obj(
        dmem,
        "DataFrame",
        "pandas",
        lambda obj: _infer_dataframe_suffix(format),
    )
    if has_dataframe:
        return dataframe_suffix

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
        lambda obj: _infer_spatialdata_suffix(format),
    )
    if has_spatialdata:
        return spatialdata_suffix
    else:
        raise NotImplementedError


def _infer_anndata_suffix(format: str | dict[str, Any] | None) -> str:
    assert not isinstance(format, dict)  # noqa: S101
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


def _infer_dataframe_suffix(format: str | dict[str, Any] | None) -> str:
    if isinstance(format, str):
        if format == ".csv":
            return ".csv"
    elif isinstance(format, dict):
        if format.get("suffix") == ".csv":
            return ".csv"
    return ".parquet"


def _infer_spatialdata_suffix(format: str | dict[str, Any] | None) -> str:
    if format is None:
        return ".zarr"
    if isinstance(format, str) and format in {"spatialdata.zarr", "zarr"}:
        return format
    raise ValueError(
        "Error when specifying SpatialData storage format, it should be"
        f" 'zarr', 'spatialdata.zarr', not '{format}'. Check 'format'"
        " or the suffix of 'key'."
    )


def write_to_disk(dmem: SupportedDataTypes, filepath: AnyPathStr, **kwargs) -> None:
    """Writes the passed in memory data to disk to a specified path."""
    if with_package_obj(
        dmem,
        "AnnData",
        "anndata",
        lambda obj: _write_anndata(obj, filepath, **kwargs),
    )[0]:
        return

    if with_package_obj(
        dmem,
        "DataFrame",
        "pandas",
        lambda obj: _write_dataframe(obj, filepath, **kwargs),
    )[0]:
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


def _write_anndata(dmem: Any, filepath: AnyPathStr, **kwargs) -> None:
    # Path(UPath(...)) properly coerces local UPaths and throws an error for cloud UPaths
    suffix = Path(filepath).suffix
    if suffix == ".h5ad":
        dmem.write_h5ad(filepath, **kwargs)
        return
    elif suffix == ".zarr":
        dmem.write_zarr(filepath, **kwargs)
        return
    else:
        raise NotImplementedError


def _write_dataframe(dmem: Any, filepath: AnyPathStr, **kwargs) -> None:
    # Path(UPath(...)) properly coerces local UPaths and throws an error for cloud UPaths
    suffix = Path(filepath).suffix
    if suffix == ".csv":
        dmem.to_csv(filepath, **kwargs)
        return
    dmem.to_parquet(filepath, **kwargs)
