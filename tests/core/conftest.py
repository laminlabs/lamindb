import shutil
from pathlib import Path
from subprocess import DEVNULL, run
from time import perf_counter

import anndata as ad
import lamindb as ln
import lamindb_setup as ln_setup
import numpy as np
import pandas as pd
import pytest

# for artifact fixtures
import yaml  # type: ignore
from lamin_utils import logger
from laminci.db import setup_local_test_postgres


def pytest_sessionstart():
    t_execute_start = perf_counter()

    ln_setup._TESTING = True
    try:
        pgurl = setup_local_test_postgres()
    except RuntimeError:
        run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
        pgurl = setup_local_test_postgres()
    ln.setup.init(
        storage="./default_storage_unit_core",
        modules="bionty",
        name="lamindb-unit-tests-core",
        db=pgurl,
    )
    ln.settings.creation.artifact_silence_missing_run_warning = True
    total_time_elapsed = perf_counter() - t_execute_start
    print(f"time to setup the instance: {total_time_elapsed:.1f}s")


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    shutil.rmtree("./default_storage_unit_core")
    ln.setup.delete("lamindb-unit-tests-core", force=True)
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602


@pytest.fixture
def ccaplog(caplog):
    """Add caplog handler to our custom logger at session start."""
    from lamin_utils._logger import logger

    logger.addHandler(caplog.handler)

    yield caplog

    logger.removeHandler(caplog.handler)


@pytest.fixture(
    scope="function",
    params=[
        # tuple of is_in_registered_storage, path, suffix, hash of test_dir
        (True, "./default_storage_unit_core/", ".csv", "iGtHiFEBV3r1_TFovdQCgw"),
        (True, "./default_storage_unit_core/", "", "iGtHiFEBV3r1_TFovdQCgw"),
        (True, "./registered_storage/", ".csv", "iGtHiFEBV3r1_TFovdQCgw"),
        (True, "./registered_storage/", "", "iGtHiFEBV3r1_TFovdQCgw"),
        (False, "./nonregistered_storage/", ".csv", "iGtHiFEBV3r1_TFovdQCgw"),
        (False, "./nonregistered_storage/", "", "iGtHiFEBV3r1_TFovdQCgw"),
    ],
)
def get_test_filepaths(request):  # -> Tuple[bool, Path, Path, Path, str]
    is_in_registered_storage: bool = request.param[0]
    root_dir: Path = Path(request.param[1])
    suffix: str = request.param[2]
    hash_test_dir: str = request.param[3]
    if is_in_registered_storage:
        # ensure that it's actually registered
        if ln.Storage.filter(root=root_dir.resolve().as_posix()).one_or_none() is None:
            ln.Storage(root=root_dir.resolve().as_posix(), type="local").save()
    else:
        assert (
            ln.Storage.filter(root=root_dir.resolve().as_posix()).one_or_none() is None
        )
    test_dirpath = root_dir / "my_dir/"
    test_dirpath.mkdir(parents=True, exist_ok=True)
    # create a first file
    test_filepath0 = test_dirpath / f"my_file{suffix}"
    test_filepath0.write_text("0")
    # create a second, duplicated file
    test_filepath1 = test_dirpath / f"my_file1{suffix}"
    test_filepath1.write_text("0")
    # create a non-duplicated file
    test_filepath2 = test_dirpath / f"my_file2{suffix}"
    test_filepath2.write_text("1")
    # return a boolean indicating whether test filepath is in default storage
    # and the test filepath
    yield (
        is_in_registered_storage,
        root_dir,
        test_dirpath,
        test_filepath0,
        suffix,
        hash_test_dir,
    )
    shutil.rmtree(test_dirpath)


@pytest.fixture(scope="session")
def example_dataframe():
    return pd.DataFrame({"feat1": [1, 2], "feat2": [3, 4]})


@pytest.fixture(scope="session")
def adata_file():
    adata = ad.AnnData(
        X=np.array([[1, 2, 3], [4, 5, 6]]),
        obs={"feat1": ["A", "B"]},
        var=pd.DataFrame(index=["MYC", "TCF7", "GATA1"]),
        obsm={"X_pca": np.array([[1, 2], [3, 4]])},
    )
    filepath = Path("adata_file.h5ad")
    adata.write(filepath)
    yield "adata_file.h5ad"
    filepath.unlink()


@pytest.fixture(scope="session")
def tsv_file():
    filepath = Path("test.tsv")
    pd.DataFrame([1, 2]).to_csv(filepath, sep="\t")
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="session")
def zip_file():
    filepath = Path("test.zip")
    pd.DataFrame([1, 2]).to_csv(filepath, sep="\t")
    yield filepath
    filepath.unlink(missing_ok=True)


@pytest.fixture(scope="session")
def yaml_file():
    filepath = Path("test.yaml")
    dct = {"a": 1, "b": 2}
    with open(filepath, "w") as f:
        yaml.dump(dct, f)
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="session")
def fcs_file():
    fcs_path = ln.examples.datasets.file_fcs_alpert19()
    yield fcs_path
    fcs_path.unlink()


@pytest.fixture(scope="session")
def mudata_file(get_small_mdata):
    filepath = Path("test.h5mu")
    get_small_mdata.write(filepath)
    yield filepath
    filepath.unlink()


@pytest.fixture(scope="session")
def spatialdata_file(get_small_sdata):
    filepath = Path("test.zarr")
    get_small_sdata.write(filepath)
    yield filepath
    shutil.rmtree(filepath)
