import os
import shutil
from pathlib import Path
from time import perf_counter

import lamindb as ln
import lamindb_setup as ln_setup
import numpy as np
import pandas as pd
import pytest
from lamin_utils import logger


def pytest_sessionstart():
    t_execute_start = perf_counter()
    ln_setup._TESTING = True
    os.environ["LAMIN_TESTING"] = "true"
    os.environ["LAMINDB_TEST_DB_VENDOR"] = "sqlite"
    print("running tests on SQLite")
    ln.setup.init(
        storage="./default_storage_tiledbsoma",
        modules="bionty",
        name="lamindb-unit-tests-tiledbsoma",
    )
    ln.settings.creation.artifact_silence_missing_run_warning = True
    # Pre-register remote roots used in tests so `ln.settings.storage = ...`
    # doesn't prompt for interactive confirmation under pytest capture.
    ln.Storage("s3://lamindb-test/tiledbsoma").save()
    total_time_elapsed = perf_counter() - t_execute_start
    print(f"time to setup the instance: {total_time_elapsed:.1f}s")


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    if Path("./default_storage_tiledbsoma").exists():
        shutil.rmtree("./default_storage_tiledbsoma")
    upath = ln_setup.core.upath.UPath("s3://lamindb-test/tiledbsoma")
    if upath.exists():
        upath.rmdir()
    ln.setup.delete("lamindb-unit-tests-tiledbsoma", force=True)
    del os.environ["LAMIN_TESTING"]


@pytest.fixture(scope="session")
def adata_file():
    import anndata as ad

    adata = ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )
    filepath = Path("adata_file.h5ad")
    adata.write(filepath)
    yield "adata_file.h5ad"
    filepath.unlink(missing_ok=True)


@pytest.fixture(scope="function")
def clean_soma_files(request):
    path = request.param if hasattr(request, "param") else "small_dataset.tiledbsoma"
    if Path(path).exists():
        shutil.rmtree(path)

    yield path

    if Path(path).exists():
        shutil.rmtree(path)


@pytest.fixture(scope="function")
def soma_experiment_file(clean_soma_files):
    import tiledbsoma.io

    adata = ln.examples.datasets.mini_immuno.get_dataset1(otype="AnnData")
    tiledbsoma.io.from_anndata("test.tiledbsoma", adata, measurement_name="RNA")
    yield "test.tiledbsoma"
    if Path("test.tiledbsoma").exists():
        shutil.rmtree("test.tiledbsoma")
