import shutil
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
