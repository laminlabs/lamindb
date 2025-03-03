import json

import psycopg2


def sign_jwt(db_url, payload: dict) -> str:
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
