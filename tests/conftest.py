import lamindb_setup
import pytest


def pytest_sessionstart(session: pytest.Session):
    if lamindb_setup._USE_DJANGO:
        lamindb_setup.init(storage="./lamindb-unit-tests")
    else:
        lamindb_setup.init(storage="./lamindb-unit-tests", schema="bionty")


def pytest_sessionfinish(session: pytest.Session):
    lamindb_setup.delete("lamindb-unit-tests")
