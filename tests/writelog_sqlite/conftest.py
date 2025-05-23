import shutil
from time import perf_counter

import lamindb_setup as ln_setup
import pytest

AUTO_CONNECT = ln_setup.settings.auto_connect
ln_setup.settings.auto_connect = False

import lamindb as ln


def pytest_sessionstart():
    t_execute_start = perf_counter()

    ln_setup._TESTING = True
    ln.setup.init(
        storage="./default_storage_writelog_core",
        modules="bionty",
        name="lamindb-writelog-tests-core",
        db=None,
    )
    ln.setup.settings.auto_connect = True
    ln.settings.creation.artifact_silence_missing_run_warning = True

    total_time_elapsed = perf_counter() - t_execute_start
    print(f"Time to setup the instance: {total_time_elapsed:.3f}s")


def pytest_sessionfinish(session: pytest.Session):
    shutil.rmtree("./default_storage_writelog_core")

    ln.setup.delete("lamindb-writelog-tests-core", force=True)
    ln.setup.settings.auto_connect = AUTO_CONNECT
