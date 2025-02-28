import lamindb as ln
from utils import _sign_jwt

pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url
token = _sign_jwt(pgurl, {"account_id": ln.setup.settings.user._uuid.hex})
# _set_token(token)


def test_fine_grained_permissions():
    assert ln.ULabel.filter().count() == 2
