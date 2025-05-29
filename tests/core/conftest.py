import shutil
from pathlib import Path
from subprocess import DEVNULL, run
from time import perf_counter

import lamindb_setup as ln_setup
import pytest
from lamin_utils import logger
from laminci.db import setup_local_test_postgres

AUTO_CONNECT = ln_setup.settings.auto_connect
ln_setup.settings.auto_connect = False

import lamindb as ln


def pytest_sessionstart():
    t_execute_start = perf_counter()

    ln_setup._TESTING = True
    pgurl = setup_local_test_postgres()
    ln.setup.init(
        storage="./default_storage_unit_core",
        modules="bionty",
        name="lamindb-unit-tests-core",
        db=pgurl,
    )
    ln.setup.settings.auto_connect = True
    ln.settings.creation.artifact_silence_missing_run_warning = True
    total_time_elapsed = perf_counter() - t_execute_start
    print(f"Time to setup the instance: {total_time_elapsed:.3f}s")


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    shutil.rmtree("./default_storage_unit_core")
    ln.setup.delete("lamindb-unit-tests-core", force=True)
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
    ln.setup.settings.auto_connect = AUTO_CONNECT


@pytest.fixture
def ccaplog(caplog):
    """Add caplog handler to our custom logger at session start."""
    from lamin_utils._logger import logger

    # Add caplog's handler to our custom logger
    logger.addHandler(caplog.handler)

    yield caplog

    # Clean up at the end of the session
    logger.removeHandler(caplog.handler)


@pytest.fixture(
    scope="module",
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
    import lamindb as ln

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
