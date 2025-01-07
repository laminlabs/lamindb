from __future__ import annotations

from pathlib import PurePosixPath
from typing import TYPE_CHECKING

from anndata import AnnData
from pandas import DataFrame

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr


def _mudata_is_installed():
    try:
        import mudata  # noqa: F401c
    except ImportError:
        return False
    return True


def infer_suffix(dmem, adata_format: str | None = None):
    """Infer LaminDB storage file suffix from a data object."""
    if isinstance(dmem, AnnData):
        if adata_format is not None:
            if adata_format not in {"h5ad", "zarr", "anndata.zarr"}:
                raise ValueError(
                    "Error when specifying AnnData storage format, it should be"
                    f" 'h5ad', 'zarr', not '{adata_format}'. Check 'format'"
                    " or the suffix of 'key'."
                )
            return "." + adata_format
        return ".h5ad"
    elif isinstance(dmem, DataFrame):
        return ".parquet"
    else:
        if _mudata_is_installed():
            from mudata import MuData

            if isinstance(dmem, MuData):
                return ".h5mu"
        raise NotImplementedError


def write_to_disk(dmem, filepath: UPathStr):
    if isinstance(dmem, AnnData):
        suffix = PurePosixPath(filepath).suffix
        if suffix == ".h5ad":
            dmem.write_h5ad(filepath)
        elif suffix == ".zarr":
            dmem.write_zarr(filepath)
        else:
            raise NotImplementedError
    elif isinstance(dmem, DataFrame):
        dmem.to_parquet(filepath)
    else:
        if _mudata_is_installed():
            from mudata import MuData

            if isinstance(dmem, MuData):
                dmem.write(filepath)
                return
        raise NotImplementedError
