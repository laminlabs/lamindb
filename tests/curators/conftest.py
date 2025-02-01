import shutil

import lamindb_setup as ln_setup
import pytest


def pytest_sessionstart():
    ln_setup.init(storage="./testdb", modules="bionty,wetlab")


def pytest_sessionfinish(session: pytest.Session):
    shutil.rmtree("./testdb")
    ln_setup.delete("testdb", force=True)
