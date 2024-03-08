import shutil
from subprocess import DEVNULL, run

import lamindb as ln
import lamindb_setup as ln_setup
import pytest
from lamin_utils import logger
from laminci.db import setup_local_test_postgres


def pytest_sessionstart():
    ln_setup._TESTING = True
    pgurl = setup_local_test_postgres()
    ln.setup.init(
        storage="./default_storage",
        schema="bionty",
        name="lamindb-unit-tests",
        db=pgurl,
    )
    ln.setup.settings.auto_connect = True


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    ln.setup.delete("lamindb-unit-tests", force=True)
    shutil.rmtree("./default_storage")
    # shutil.rmtree("./outside_storage")
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)
    ln.setup.settings.auto_connect = False
