import shutil
from subprocess import DEVNULL, run
from time import perf_counter

import lamindb_setup as ln_setup
import pytest
from lamin_utils import logger


def pytest_sessionstart():
    t_execute_start = perf_counter()

    ln_setup.settings.auto_connect = True
    # these are called in separate scripts because can't change connection
    # within the same python process due to django
    # init instance and setup RLS
    run(  # noqa: S602
        "python ./tests/permissions/scripts/setup_instance.py",
        shell=True,
        capture_output=False,
    )
    # populate permissions and models via the admin connection
    run(  # noqa: S602
        "python ./tests/permissions/scripts/setup_access.py",
        shell=True,
        capture_output=False,
    )

    total_time_elapsed = perf_counter() - t_execute_start
    print(f"Time to setup the instance: {total_time_elapsed:.3f}s")


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    shutil.rmtree("./default_storage_permissions")
    ln_setup.delete("lamindb-test-permissions", force=True)
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
