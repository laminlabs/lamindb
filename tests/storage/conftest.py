import shutil
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
    ln.settings.storage = "./default_storage_unit_storage"


def pytest_sessionstart():
    t_execute_start = perf_counter()

    ln_setup._TESTING = True
    pgurl = setup_local_test_postgres()
    try:
        create_test_instance(pgurl)
    except Exception as e:
        print("failed to create test instance:", e)
        delete_test_instance()
        print("retrying after deleting the instance")
        create_test_instance(pgurl)
    total_time_elapsed = perf_counter() - t_execute_start
    print(f"time to setup the instance: {total_time_elapsed:.1f}s")


def delete_test_instance():
    logger.set_verbosity(1)
    shutil.rmtree("./default_storage_unit_storage")
    # handle below better in the future
    if ln.UPath("s3://lamindb-test/storage/.lamindb").exists():
        ln.UPath("s3://lamindb-test/storage/.lamindb").rmdir()
    another_storage = ln.UPath("s3://lamindb-ci/lamindb-unit-tests-cloud/.lamindb")
    if another_storage.exists():
        another_storage.rmdir()
    ln.setup.delete("lamindb-unit-tests-storage", force=True)
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
    ln.setup.settings.auto_connect = AUTO_CONNECT


def pytest_sessionfinish(session: pytest.Session):
    delete_test_instance()


@pytest.fixture
def ccaplog(caplog):
    """Add caplog handler to our custom logger at session start."""
    from lamin_utils._logger import logger

    # Add caplog's handler to our custom logger
    logger.addHandler(caplog.handler)

    yield caplog

    # Clean up at the end of the session
    logger.removeHandler(caplog.handler)
