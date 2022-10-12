from typing import Optional

from anndata import AnnData
from pandas import DataFrame


def infer_suffix(dmem, adata_format: Optional[str] = None):
    """Infer LaminDB storage file suffix from a data object."""
    if isinstance(dmem, AnnData):
        if adata_format is not None:
            if adata_format not in ("h5ad", "zarr"):
                raise ValueError
            return "." + adata_format
        return ".h5ad"
    elif isinstance(dmem, DataFrame):
        return ".feather"
    else:
        raise NotImplementedError


def write_to_file(dmem, filepath: str):
    if isinstance(dmem, AnnData):
        dmem.write(filepath)
    elif isinstance(dmem, DataFrame):
        try:
            dmem.to_feather(filepath)
        except ValueError:
            dmem.reset_index().to_feather(filepath)
    else:
        raise NotImplementedError
