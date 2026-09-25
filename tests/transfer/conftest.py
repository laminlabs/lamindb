import os
import shutil
import time
from subprocess import DEVNULL, run

import lamindb as ln
import pytest
from laminci.db import setup_local_test_postgres

_POSTGRES_CONTAINER = "pgtest"


def _close_all_connections():
    from django.db import connections

    connections.close_all()


def _is_postgresql() -> bool:
    return os.getenv("LAMINDB_TEST_DB_VENDOR") == "postgresql"


def _postgres_url(database: str) -> str:
    return f"postgresql://postgres:pwd@0.0.0.0:5432/{database}"


def _init_instance(name: str) -> None:
    kwargs = {"storage": f"./{name}", "name": name}
    if _is_postgresql():
        kwargs["db"] = _postgres_url(name)
    ln.setup.init(**kwargs)


def _wait_for_postgres() -> None:
    for _ in range(30):
        process = run(
            f"docker exec {_POSTGRES_CONTAINER} pg_isready -U postgres",
            shell=True,
            stdout=DEVNULL,
            stderr=DEVNULL,
        )
        if process.returncode == 0:
            return
        time.sleep(1)
    raise RuntimeError("Postgres test container did not become ready")


def _create_database(name: str) -> None:
    run(
        f'docker exec {_POSTGRES_CONTAINER} psql -U postgres -c "CREATE DATABASE {name}"',
        shell=True,
        check=True,
    )


def pytest_sessionstart():
    if _is_postgresql():
        print("running transfer tests on PostgreSQL")
        try:
            setup_local_test_postgres()
        except RuntimeError:
            run(
                f"docker stop {_POSTGRES_CONTAINER} && docker rm {_POSTGRES_CONTAINER}",
                shell=True,
                stdout=DEVNULL,
            )
            setup_local_test_postgres()
        _wait_for_postgres()
        _create_database("testdb1")
        _create_database("testdb2")
    else:
        os.environ["LAMINDB_TEST_DB_VENDOR"] = "sqlite"
        print("running transfer tests on SQLite")


def pytest_sessionfinish(session: pytest.Session):
    if _is_postgresql():
        run(
            f"docker stop {_POSTGRES_CONTAINER} && docker rm {_POSTGRES_CONTAINER}",
            shell=True,
            stdout=DEVNULL,
        )


@pytest.fixture(scope="session", autouse=True)
def setup_testdb1():
    _init_instance("testdb1")
    yield
    _close_all_connections()
    shutil.rmtree("./testdb1")
    ln.setup.delete("testdb1", force=True)


@pytest.fixture(scope="session", autouse=True)
def setup_testdb2():
    _init_instance("testdb2")
    yield
    _close_all_connections()
    shutil.rmtree("./testdb2")
    ln.setup.delete("testdb2", force=True)


@pytest.fixture(scope="session", autouse=True)
def setup_test_in_db1_then_connect_to_db2(setup_testdb1, setup_testdb2):
    ln.connect("testdb1")
    ln.Artifact("README.md", key="README.md").save()
    ln.connect("testdb2")
