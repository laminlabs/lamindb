import shutil
from pathlib import Path

import anndata as ad
import lamindb as ln
import numpy as np
import pandas as pd
import pytest
from lamindb.core.storage._zarr import identify_zarr_type
from lamindb_setup.core.upath import (
    CloudPath,
)


@pytest.fixture(scope="session")
def get_small_adata():
    return ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )


def test_zarr_upload_cache(get_small_adata):
    previous_storage = ln.setup.settings.storage.root_as_str
    ln.settings.storage = "s3://lamindb-test/core"

    zarr_path = Path("./test_adata.zarr")
    get_small_adata.write_zarr(zarr_path)

    artifact = ln.Artifact(zarr_path, key="test_adata.zarr")
    assert artifact.otype == "AnnData"
    assert artifact.n_files >= 1
    artifact.save()

    assert isinstance(artifact.path, CloudPath)
    assert artifact.path.exists()
    assert identify_zarr_type(artifact.path) == "anndata"

    shutil.rmtree(artifact.cache())

    cache_path = artifact._cache_path
    assert isinstance(artifact.load(), ad.AnnData)
    assert cache_path.is_dir()

    shutil.rmtree(cache_path)
    assert not cache_path.exists()
    artifact.cache()
    assert cache_path.is_dir()

    artifact.delete(permanent=True, storage=True)
    shutil.rmtree(zarr_path)

    # test zarr from memory
    artifact = ln.Artifact(get_small_adata, key="test_adata.anndata.zarr")
    assert artifact._local_filepath.is_dir()
    assert artifact.otype == "AnnData"
    assert artifact.suffix == ".anndata.zarr"
    assert artifact.n_files >= 1

    artifact.save()
    assert isinstance(artifact.path, CloudPath)
    assert artifact.path.exists()
    cache_path = artifact._cache_path
    assert cache_path.is_dir()

    shutil.rmtree(cache_path)
    assert not cache_path.exists()

    artifact._memory_rep = None

    assert isinstance(artifact.load(), ad.AnnData)
    assert cache_path.is_dir()

    artifact.delete(permanent=True, storage=True)

    ln.settings.storage = previous_storage
