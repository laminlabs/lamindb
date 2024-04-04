from __future__ import annotations

from typing import TYPE_CHECKING

from anndata import AnnData
from pandas import DataFrame

if TYPE_CHECKING:
    from lamindb_setup.core.types import UPathStr


def infer_suffix(dmem, adata_format: str | None = None):
    """Infer LaminDB storage file suffix from a data object."""
    if isinstance(dmem, AnnData):
        if adata_format is not None:
            # below should be zrad, not zarr
            if adata_format not in ("h5ad", "zarr", "zrad"):
                raise ValueError(
                    "Error when specifying AnnData storage format, it should be"
                    f" 'h5ad', 'zarr' or 'zrad', not '{adata_format}'. Check 'format'"
                    " or the suffix of 'key'."
                )
            return "." + adata_format
        return ".h5ad"
    elif isinstance(dmem, DataFrame):
        return ".parquet"
    else:
        raise NotImplementedError


def write_to_file(dmem, filepath: UPathStr):
    if isinstance(dmem, AnnData):
        dmem.write(filepath)
    elif isinstance(dmem, DataFrame):
        dmem.to_parquet(filepath)
    else:
        raise NotImplementedError
