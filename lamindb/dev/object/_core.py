from pathlib import Path
from typing import Optional, Union

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
        return ".parquet"
    else:
        raise NotImplementedError


def write_to_file(dmem, filepath: Union[str, Path]):
    if isinstance(dmem, AnnData):
        dmem.write(filepath)
    elif isinstance(dmem, DataFrame):
        try:
            dmem.to_parquet(filepath)
        except ValueError:
            dmem.reset_index().to_parquet(filepath)
    else:
        raise NotImplementedError
