import lamindb as ln  # noqa
import hubmodule.models as hm
from ..utils import _create_jwt_user

full_access = ln.models.Space(name="full access", uid="00000001").save()
select_access = ln.models.Space(name="full access", uid="00000002").save()
no_access = ln.models.Space(name="no access", uid="00000003").save()

account = hm.Account(id=ln.setup.settings.user._uuid.hex).save()

# no access space
ulabel = ln.ULabel(name="no_access_ulabel")
ulabel.space = no_access
ulabel.save()
# setup full access space
hm.AccessSpace(account=account, space=full_access, operation="SELECT").save()
hm.AccessSpace(account=account, space=full_access, operation="INSERT").save()
hm.AccessSpace(account=account, space=full_access, operation="UPDATE").save()
hm.AccessSpace(account=account, space=full_access, operation="DELETE").save()

ulabel = ln.ULabel(name="full_access_ulabel")
ulabel.space = full_access
ulabel.save()
# setup select access space
hm.AccessSpace(account=account, space=select_access, operation="SELECT").save()

ulabel = ln.ULabel(name="select_ulabel")
ulabel.space = select_access
ulabel.save()

# create a db connection url that works with RLS
pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url
jwt_db_url = _create_jwt_user(pgurl)  #
ln.setup.settings.instance._db = jwt_db_url
ln.setup.settings.instance._persist()
