import os
import shutil
from pathlib import Path

import lamindb as ln
import lamindb_setup as ln_setup
import pytest
from lamindb_setup.core._settings_store import local_current_instance_file

INSTANCE_NAME = "lamindb-unit-tests-minimal"
STORAGE_ROOT = "./default_storage_minimal"
# Directories created by blob constructor fixtures, outside the instance storage.
FIXTURE_STORAGE_ROOTS = ("./registered_storage", "./nonregistered_storage")


def pytest_sessionstart():
    ln_setup._TESTING = True
    os.environ["LAMIN_TESTING"] = "true"
    os.environ["LAMINDB_TEST_DB_VENDOR"] = "sqlite"
    # Avoid auto-connecting to a stale local instance marker in the repo root.
    local_current_instance_file(Path.cwd()).unlink(missing_ok=True)
    storage = Path(STORAGE_ROOT)
    if storage.exists():
        shutil.rmtree(storage)
    ln.setup.init(storage=STORAGE_ROOT, name=INSTANCE_NAME)
    ln.settings.creation.artifact_silence_missing_run_warning = True


def pytest_sessionfinish():
    ln_setup.settings.dev_dir = None
    shutil.rmtree(STORAGE_ROOT, ignore_errors=True)
    for root in FIXTURE_STORAGE_ROOTS:
        shutil.rmtree(root, ignore_errors=True)
    ln.setup.delete(INSTANCE_NAME, force=True)


@pytest.fixture(
    scope="function",
    params=[
        # tuple of is_in_registered_storage, path, suffix, hash of test_dir
        (True, f"{STORAGE_ROOT}/", ".csv", "iGtHiFEBV3r1_TFovdQCgw"),
        (True, f"{STORAGE_ROOT}/", "", "iGtHiFEBV3r1_TFovdQCgw"),
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
