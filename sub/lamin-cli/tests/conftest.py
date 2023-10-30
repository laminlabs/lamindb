import shutil

import lamindb_setup
import pytest
from lamin_utils import logger


def pytest_sessionstart(session: pytest.Session):
    lamindb_setup.init(
        storage="./default_storage",
        name="lamindb-setup-notebook-tests",
    )
    lamindb_setup._TESTING = True


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    lamindb_setup.delete("lamindb-setup-notebook-tests", force=True)
    shutil.rmtree("./default_storage")
