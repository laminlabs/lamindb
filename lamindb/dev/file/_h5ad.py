import anndata
import fsspec
from anndata import AnnData, read_h5ad
from lndb_setup import settings


def h5ad_to_anndata(filekey) -> AnnData:
    """h5ad → AnnData."""
    return anndata.read(settings.instance.storage.local_filepath(filekey))


def read_adata_h5ad(filepath, **kwargs) -> AnnData:
    with fsspec.open(filepath, mode="rb") as file:
        adata = read_h5ad(file, backed=False, **kwargs)
        return adata
