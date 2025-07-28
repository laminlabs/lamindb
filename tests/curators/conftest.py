import shutil
from time import perf_counter

import lamindb_setup as ln_setup
import pytest


def pytest_sessionstart():
    t_execute_start = perf_counter()
    ln_setup.init(storage="./testdb", modules="bionty,wetlab")
    total_time_elapsed = perf_counter() - t_execute_start
    print(f"time to setup the instance: {total_time_elapsed:.1f}s")


def pytest_sessionfinish(session: pytest.Session):
    shutil.rmtree("./testdb")
    ln_setup.delete("testdb", force=True)


@pytest.fixture
def ccaplog(caplog):
    """Add caplog handler to our custom logger at session start."""
    from lamin_utils._logger import logger

    # Add caplog's handler to our custom logger
    logger.addHandler(caplog.handler)

    yield caplog

    # Clean up at the end of the session
    logger.removeHandler(caplog.handler)
