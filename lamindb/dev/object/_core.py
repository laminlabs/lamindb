from anndata import AnnData
from pandas import DataFrame


def infer_file_suffix(dmem):
    """Infer LaminDB storage file suffix from a data object."""
    if isinstance(dmem, AnnData):
        return ".h5ad"
    elif isinstance(dmem, DataFrame):
        return ".feather"
    else:
        raise NotImplementedError
