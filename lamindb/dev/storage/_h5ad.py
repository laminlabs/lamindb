import anndata
from anndata import AnnData
from lamindb_setup import settings
from lamindb_setup.dev.upath import infer_filesystem


def h5ad_to_anndata(filekey) -> AnnData:
    """h5ad → AnnData."""
    return anndata.read(settings.instance.storage.local_filepath(filekey))


def read_adata_h5ad(filepath, **kwargs) -> AnnData:
    fs, filepath = infer_filesystem(filepath)

    with fs.open(filepath, mode="rb") as file:
        adata = anndata.read_h5ad(file, backed=False, **kwargs)
        return adata
