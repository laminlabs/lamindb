import shutil

import lamindb_setup
import pytest


def pytest_sessionstart(session: pytest.Session):
    lamindb_setup.init(
        storage="./default_storage", schema="bionty", name="lamindb-unit-tests"
    )


def pytest_sessionfinish(session: pytest.Session):
    lamindb_setup.delete("lamindb-unit-tests")
    shutil.rmtree("./default_storage")
