import json

import psycopg2
from laminhub_rest.core.db import DbRoleHandler


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


def _create_jwt_user(dsn_admin: str):
    db_role_handler = DbRoleHandler(dsn_admin)
    jwt_role_name = "permissions_jwt"
    jwt_db_url = db_role_handler.create(
        jwt_role_name, expires_in=None, alter_if_exists=True
    )
    db_role_handler.permission.grant_write_jwt(jwt_role_name)
    return jwt_db_url
