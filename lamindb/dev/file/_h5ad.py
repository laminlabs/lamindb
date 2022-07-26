import anndata
from anndata import AnnData
from lndb_setup import settings


def h5ad_to_anndata(filekey) -> AnnData:
    """h5ad → AnnData."""
    return anndata.read(settings.storage.local_filepath(filekey))
