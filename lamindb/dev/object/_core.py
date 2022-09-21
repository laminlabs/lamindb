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
