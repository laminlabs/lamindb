import json
import shutil
from subprocess import DEVNULL, run
from time import perf_counter

import lamindb_setup as ln_setup
import psycopg2
import pytest
from lamin_utils import logger
from laminci.db import setup_local_test_postgres
from laminhub_rest.core.db import DbRoleHandler
from laminhub_rest.hubmodule.hubmodule._setup import _install_db_module


def _create_jwt_user(dsn_admin: str):
    db_role_handler = DbRoleHandler(dsn_admin)
    jwt_role_name = "permissions_jwt"
    jwt_db_url = db_role_handler.create(
        jwt_role_name, expires_in=None, alter_if_exists=True
    )
    db_role_handler.permission.grant_write_jwt(jwt_role_name)
    return jwt_db_url


def _sign_jwt(db_url, payload: dict) -> str:
    with psycopg2.connect(db_url) as conn, conn.cursor() as cur:
        cur.execute(
            """
                SELECT sign(
                    %s::json,
                    (SELECT security.get_secret('jwt_secret')),
                    %s
                )
                """,
            (json.dumps(payload), "HS256"),
        )
        token = cur.fetchone()[0]
        if not token:
            msg = "Failed to generate JWT"
            raise ValueError(msg)
        return token


def pytest_sessionstart():
    t_execute_start = perf_counter()

    pgurl = setup_local_test_postgres()
    ln_setup.init(
        storage="./default_storage_permissions",
        name="lamindb-test-permissions",
        db=pgurl,
    )
    _install_db_module(pgurl)
    jwturl = _create_jwt_user(pgurl)
    print(jwturl)

    total_time_elapsed = perf_counter() - t_execute_start
    print(f"Time to setup the instance: {total_time_elapsed:.3f}s")


def pytest_sessionfinish(session: pytest.Session):
    logger.set_verbosity(1)
    shutil.rmtree("./default_storage_permissions")
    ln_setup.delete("lamindb-test-permissions", force=True)
    run("docker stop pgtest && docker rm pgtest", shell=True, stdout=DEVNULL)  # noqa: S602
