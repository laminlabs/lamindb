from django.db import connection


def test_connection_is_sqlite():
    assert connection.vendor == "sqlite"
