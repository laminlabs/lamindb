import shutil
from pathlib import Path
from typing import Generator

import lamindb as ln
import numpy as np
import pandas as pd
import pytest
from scipy.sparse import csr_matrix


@pytest.fixture(scope="session")
def get_small_adata():
    # shouldn't need anndata installed to run tests
    import anndata as ad

    return ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )


@pytest.fixture(scope="session")
def get_small_mdata():
    # shouldn't need mudata installed to run tests
    import anndata as ad
    import mudata as md

    adata1 = ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )

    adata2 = ad.AnnData(
        X=np.array([[7, 8], [9, 10]]),
        obs={"feat2": ["C", "D"]},
        var=pd.DataFrame(index=["FOXP3", "CD8A"]),
        obsm={"X_umap": np.array([[5, 6], [7, 8]])},
    )

    return md.MuData({"rna": adata1, "protein": adata2})


@pytest.fixture(scope="session")
def get_small_sdata():
    # shouldn't need spatialdata installed to run tests
    import anndata as ad
    import spatialdata as sd

    adata = ad.AnnData(
        X=csr_matrix(np.array([[0.1, 0.2], [0.3, 0.4]])),
        obs=pd.DataFrame(index=["cell1", "cell2"]),
        var=pd.DataFrame(index=["gene1", "gene2"]),
    )

    {
        "region1": np.array([[[0, 0], [0, 1], [1, 1], [1, 0]]]),
        "region2": np.array([[[2, 2], [2, 3], [3, 3], [3, 2]]]),
    }

    sdata_obj = sd.SpatialData(
        tables={"gene_expression": adata},
    )

    return sdata_obj


@pytest.fixture(scope="session")
def get_small_soma_experiment():
    # shouldn't need tilebdomsa installed to run tests
    import tiledbsoma
    import tiledbsoma.io

    adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
    tiledbsoma.io.from_anndata("test.tiledbsoma", adata, measurement_name="RNA")

    exp = tiledbsoma.Experiment.open("test.tiledbsoma")
    yield exp

    shutil.rmtree("test.tiledbsoma")


@pytest.fixture(scope="session")
def get_mini_csv() -> Generator[Path, None, None]:
    csv_path = ln.examples.datasets.file_mini_csv()
    yield csv_path

    Path("mini.csv").unlink(missing_ok=True)
