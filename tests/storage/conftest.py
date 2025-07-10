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


def create_test_instance(pgurl: str):
    ln.setup.init(
        storage="./default_storage_unit_storage",
        modules="bionty",
        name="lamindb-unit-tests-storage",
        db=pgurl,
    )
    ln.setup.register()  # temporarily
    ln.setup.settings.auto_connect = True
    ln.settings.creation.artifact_silence_missing_run_warning = True
    ln.settings.storage = (
        "s3://lamindb-ci/test-data"  # register as valid storage location
    )
    ln.settings.storage = "./default_storage_unit_storage", "testuser1-laptop"


def pytest_sessionstart():
    t_execute_start = perf_counter()

    ln_setup._TESTING = True
    try:
        pgurl = setup_local_test_postgres()
    except RuntimeError:
        run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
        pgurl = setup_local_test_postgres()
    try:
        create_test_instance(pgurl)
    except Exception as e:
        print("failed to create test instance:", e)
        print("deleting the instance")
        delete_test_instance()
        # below currently fails because cannot create two instances in the same session
        # create_test_instance(pgurl)
        print("now rerun")
        quit()
    total_time_elapsed = perf_counter() - t_execute_start
    print(f"time to setup the instance: {total_time_elapsed:.1f}s")
    assert ln.Storage.filter(root="s3://lamindb-ci/test-data").one_or_none() is not None


def delete_test_instance():
    logger.set_verbosity(1)
    if Path("./default_storage_unit_storage").exists():
        shutil.rmtree("./default_storage_unit_storage")
    # handle below better in the future
    for path in (
        "s3://lamindb-test/storage/.lamindb",
        "s3://lamindb-test/core/.lamindb",
        "s3://lamindb-ci/lamindb-unit-tests-cloud/.lamindb",
        "s3://lamindb-ci/test-settings-switch-storage/.lamindb",
    ):
        upath = ln.UPath(path)
        if upath.exists():
            upath.rmdir()
    ln.setup.delete("lamindb-unit-tests-storage", force=True)


def pytest_sessionfinish(session: pytest.Session):
    delete_test_instance()
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
