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
    ln.setup.register()  # temporarily
    ln.setup.settings.auto_connect = True
    ln.settings.creation.artifact_silence_missing_run_warning = True

    total_time_elapsed = perf_counter() - t_execute_start
    print(f"Time to setup the instance: {total_time_elapsed:.3f}s")


def pytest_sessionfinish(session: pytest.Session):
    try:
        logger.set_verbosity(1)
        shutil.rmtree("./default_storage_unit_core")

        ln.setup.delete("lamindb-unit-tests-core", force=True)
        ln.setup.settings.auto_connect = AUTO_CONNECT
    finally:
        run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
