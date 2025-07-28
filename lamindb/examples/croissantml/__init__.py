from pathlib import Path


def mini_immuno() -> tuple[Path, Path]:
    """Return paths to the mini immuno dataset and its metadata as a CroissantML file."""
    from ...core.datasets.mini_immuno import get_dataset1

    adata = get_dataset1(otype="AnnData")
    dataset_path = Path("mini_immuno.anndata.zarr")
    adata.write_zarr(dataset_path)
    return dataset_path, Path(
        __file__
    ).parent / "mini_immuno.anndata.zarr_metadata.json"
