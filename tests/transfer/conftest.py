import os
import shutil
from subprocess import DEVNULL, run

import lamindb as ln
import pytest
from laminci.db import setup_local_test_postgres

pgurls: dict[str, str] = {}


def _close_all_connections():
    from django.db import connections

    connections.close_all()


def pytest_sessionstart():
    is_postgresql = os.getenv("LAMINDB_TEST_DB_VENDOR") == "postgresql"
    if is_postgresql:
        print("running transfer tests on PostgreSQL")
        try:
            pgurls.update(
                setup_local_test_postgres(
                    databases=["testdb1", "testdb2", "testdbbionty"]
                )
            )
        except RuntimeError:
            run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
            pgurls.update(
                setup_local_test_postgres(
                    databases=["testdb1", "testdb2", "testdbbionty"]
                )
            )
    else:
        os.environ["LAMINDB_TEST_DB_VENDOR"] = "sqlite"
        print("running transfer tests on SQLite")


def pytest_sessionfinish(session: pytest.Session):
    if os.getenv("LAMINDB_TEST_DB_VENDOR") != "sqlite":
        run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602


@pytest.fixture(scope="session", autouse=True)
def setup_testdb1():
    if os.getenv("LAMINDB_TEST_DB_VENDOR") == "postgresql":
        ln.setup.init(
            storage="./testdb1",
            name="testdb1",
            db=pgurls["testdb1"],
        )
    else:
        ln.setup.init(storage="./testdb1", name="testdb1")
    yield
    _close_all_connections()
    shutil.rmtree("./testdb1")
    ln.setup.delete("testdb1", force=True)


@pytest.fixture(scope="session", autouse=True)
def setup_testdb2():
    if os.getenv("LAMINDB_TEST_DB_VENDOR") == "postgresql":
        ln.setup.init(
            storage="./testdb2",
            name="testdb2",
            db=pgurls["testdb2"],
        )
    else:
        ln.setup.init(storage="./testdb2", name="testdb2")
    yield
    _close_all_connections()
    shutil.rmtree("./testdb2")
    ln.setup.delete("testdb2", force=True)


@pytest.fixture(scope="session", autouse=True)
def setup_test_in_db1_then_connect_to_db2(setup_testdb1, setup_testdb2):
    ln.connect("testdb1")
    ln.Artifact("README.md", key="README.md").save()
    ln.connect("testdb2")


@pytest.fixture(scope="session")
def bionty_instance(setup_testdb2):
    """Target instance with the bionty module. testdb2 stays without it."""
    name = "testdbbionty"
    kwargs = {"storage": f"./{name}", "name": name, "modules": "bionty"}
    if os.getenv("LAMINDB_TEST_DB_VENDOR") == "postgresql":
        kwargs["db"] = pgurls[name]
    ln.setup.init(**kwargs)
    ln.connect("testdb2")
    yield name
    ln.connect("testdb2")
    _close_all_connections()
    shutil.rmtree(f"./{name}")
    ln.setup.delete(name, force=True)


@pytest.fixture
def connected_bionty(bionty_instance):
    ln.connect(bionty_instance)
    yield
    ln.connect("testdb2")
