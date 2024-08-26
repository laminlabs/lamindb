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
        storage="./default_storage_unit_storage",
        name="lamindb-unit-tests-storage",
        db=pgurl,
    )
    ln.setup.register()  # temporarily
    ln.setup.settings.auto_connect = True
    ln.settings.creation.artifact_silence_missing_run_warning = True


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    shutil.rmtree("./default_storage_unit_storage")
    # handle below better in the future
    if ln.UPath("s3://lamindb-test/storage/.lamindb").exists():
        ln.UPath("s3://lamindb-test/storage/.lamindb").rmdir()
    ln.setup.delete("lamindb-unit-tests-storage", force=True)
    # shutil.rmtree("./outside_storage")
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
    ln.setup.settings.auto_connect = AUTO_CONNECT
