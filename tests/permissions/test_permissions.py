import lamindb as ln
from utils import _create_jwt_user

pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url
jwt_db_url = _create_jwt_user(pgurl)  # create a db connection url that works with RLS
ln.connect(_db=jwt_db_url)


def test_fine_grained_permissions():
    # that is just to check that everything was setup properly
    # this doesn't check the real permissions functionality
    assert ln.setup.settings.instance.modules == {"hubmodule"}

    assert ln.ULabel.filter().count() == 2
