import shutil
from pathlib import Path
from subprocess import DEVNULL, run

import lamindb_setup
import pytest
from lamin_utils import logger
from laminci.db import setup_local_test_postgres


def pytest_sessionstart(session: pytest.Session):
    pgurl = setup_local_test_postgres()
    lamindb_setup.init(
        storage="./default_storage",
        schema="bionty",
        name="lamindb-unit-tests",
        db=pgurl,
    )
    # we're setting this to true prior to importing lamindb!
    lamindb_setup._TESTING = True


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    lamindb_setup.delete("lamindb-unit-tests", force=True)
    shutil.rmtree("./default_storage")
    # shutil.rmtree("./outside_storage")
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)


@pytest.fixture(
    scope="module",
    params=[
        # tuple of isin_existing_storage, path, suffix
        (True, "./default_storage/", ".csv"),
        (True, "./default_storage/", ""),
        (True, "./registered_storage/", ".csv"),
        (True, "./registered_storage/", ""),
        (False, "./nonregistered_storage/", ".csv"),
        (False, "./nonregistered_storage/", ""),
    ],
)
def get_test_filepaths(request):  # -> Tuple[bool, Path, Path, Path, str]
    import lamindb as ln

    isin_existing_storage: bool = request.param[0]
    root_dir: Path = Path(request.param[1])
    suffix: str = request.param[2]
    if isin_existing_storage:
        # ensure that it's actually registered
        if ln.Storage.filter(root=root_dir.resolve().as_posix()).one_or_none() is None:
            ln.Storage(root=root_dir.resolve().as_posix(), type="local").save()
    test_dir = root_dir / "my_dir/"
    test_dir.mkdir(parents=True)
    test_filepath = test_dir / f"my_file{suffix}"
    test_filepath.write_text(str(test_filepath))
    # create a duplicated file
    test_filepath1 = test_dir / f"my_file1{suffix}"
    test_filepath1.write_text(str(test_filepath))
    # create a non-duplicated file
    test_filepath2 = test_dir / f"my_file2{suffix}"
    test_filepath2.write_text(str(test_filepath2))
    # return a boolean indicating whether test filepath is in default storage
    # and the test filepath
    yield (isin_existing_storage, root_dir, test_dir, test_filepath, suffix)
    shutil.rmtree(test_dir)
