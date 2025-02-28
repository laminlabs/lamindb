import lamindb as ln
from jwt_utils import set_jwt, sign_jwt

pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url
token = sign_jwt(pgurl, {"account_id": ln.setup.settings.user._uuid.hex})
set_jwt(token)


def test_fine_grained_permissions():
    assert ln.ULabel.filter().count() == 2
    # should delete
    ln.ULabel.get(name="full_access_ulabel").delete()
    assert ln.ULabel.filter().count() == 1
    # should not delete, doesn't error for some reason
    ln.ULabel.get(name="select_ulabel").delete()
    assert ln.ULabel.filter().count() == 1
