from pathlib import Path


def mini_immuno():
    from ...core.datasets.mini_immuno import get_dataset1

    adata = get_dataset1(otype="AnnData")
    adata.write("mini_immuno.anndata.zarr")
    return Path(__file__).parent / "mini_immuno.anndata.zarr_metadata.json"
