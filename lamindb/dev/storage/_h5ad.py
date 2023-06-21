import anndata
from anndata import AnnData
from lamindb_setup.dev.upath import infer_filesystem


def read_adata_h5ad(filepath, **kwargs) -> AnnData:
    fs, filepath = infer_filesystem(filepath)

    with fs.open(filepath, mode="rb") as file:
        adata = anndata.read_h5ad(file, backed=False, **kwargs)
        return adata
