import lamindb as ln
from utils import _create_jwt_user

pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"
jwt_db_url = _create_jwt_user(pgurl)


def test_fine_grained_permissions():
    ln.connect("lamindb-test-permissions", _db=jwt_db_url)

    assert ln.setup.settings.instance.modules == {"hubmodule"}
