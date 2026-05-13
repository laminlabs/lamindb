import shutil

import lamindb as ln
import pytest


@pytest.fixture(scope="session", autouse=True)
def setup_testdb1():
    ln.setup.init(storage="./testdb1")
    yield
    shutil.rmtree("./testdb1")
    ln.setup.delete("testdb1", force=True)


@pytest.fixture(scope="session", autouse=True)
def setup_testdb2():
    ln.setup.init(storage="./testdb2")
    yield
    shutil.rmtree("./testdb2")
    ln.setup.delete("testdb2", force=True)
