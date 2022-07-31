from pathlib import Path

from anndata import AnnData
from pandas import DataFrame

from ..object import anndata_to_h5ad


def infer_file_suffix(dmem):
    """Infer LaminDB storage file suffix from a data object."""
    if isinstance(dmem, AnnData):
        return ".h5ad"
    elif isinstance(dmem, DataFrame):
        return ".feather"
    else:
        raise NotImplementedError


def write_to_file(dmem, filekey: str) -> Path:
    if isinstance(dmem, AnnData):
        return anndata_to_h5ad(adata=dmem, filekey=filekey)
    else:
        raise NotImplementedError
