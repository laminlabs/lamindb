import shutil
from subprocess import run

import lamindb_setup
import pytest
from laminci.db import setup_local_test_postgres


def pytest_sessionstart(session: pytest.Session):
    pgurl = setup_local_test_postgres()
    lamindb_setup.init(
        storage="./default_storage",
        schema="bionty",
        name="lamindb-unit-tests",
        db=pgurl,
    )


def pytest_sessionfinish(session: pytest.Session):
    lamindb_setup.delete("lamindb-unit-tests")
    shutil.rmtree("./default_storage")
    run("docker stop pgtest && docker rm pgtest", shell=True)
