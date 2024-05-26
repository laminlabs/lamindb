import shutil
from pathlib import Path
from subprocess import DEVNULL, run

import lamindb as ln
import lamindb_setup as ln_setup
import pytest
from lamin_utils import logger
from laminci.db import setup_local_test_postgres

AUTO_CONNECT = ln.setup.settings.auto_connect


def pytest_sessionstart():
    ln_setup._TESTING = True
    pgurl = setup_local_test_postgres()
    ln.setup.init(
        storage="./default_storage",
        schema="bionty",
        name="lamindb-unit-tests",
        db=pgurl,
    )
    # ln.setup.register()  # temporarily
    ln.setup.settings.auto_connect = True
    ln.settings.silence_file_run_transform_warning = True


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    shutil.rmtree("./default_storage")
    # handle below better in the future
    if ln.UPath("s3://lamindb-test/.lamindb").exists():
        ln.UPath("s3://lamindb-test/.lamindb").rmdir()
    ln.setup.delete("lamindb-unit-tests", force=True)
    # shutil.rmtree("./outside_storage")
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)
    ln.setup.settings.auto_connect = AUTO_CONNECT


@pytest.fixture(
    scope="module",
    params=[
        # tuple of is_in_registered_storage, path, suffix, hash of test_dir
        (True, "./default_storage/", ".csv", "iGtHiFEBV3r1_TFovdQCgw"),
        (True, "./default_storage/", "", "iGtHiFEBV3r1_TFovdQCgw"),
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
    test_dirpath.mkdir(parents=True)
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
