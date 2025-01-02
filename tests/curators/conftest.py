import lamindb_setup as ln_setup
import pytest


def pytest_sessionstart():
    ln_setup.init(storage="./testdb", schema="bionty")


def pytest_sessionfinish(session: pytest.Session):
    ln_setup.delete("testdb", force=True)
