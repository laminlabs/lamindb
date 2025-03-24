import lamindb as ln  # noqa
import hubmodule.models as hm
from laminhub_rest.core.db import DbRoleHandler

full_access = ln.models.Space(name="full access", uid="00000001").save()  # type: ignore
select_access = ln.models.Space(name="select access", uid="00000002").save()  # type: ignore
no_access = ln.models.Space(name="no access", uid="00000003").save()  # type: ignore

account = hm.Account(id=ln.setup.settings.user._uuid.hex).save()

# space All
ulabel = ln.ULabel(name="space_all_ulabel").save()

# no access space
ulabel = ln.ULabel(name="no_access_ulabel")
ulabel.space = no_access
ulabel.save()

project = ln.Project(name="No_access_project")  # type: ignore
project.space = no_access
project.save()

# setup write access space
hm.AccessSpace(account=account, space=full_access, role="write").save()

ulabel = ln.ULabel(name="full_access_ulabel")
ulabel.space = full_access
ulabel.save()
# setup read access space
hm.AccessSpace(account=account, space=select_access, role="read").save()

ulabel = ln.ULabel(name="select_ulabel")
ulabel.space = select_access
ulabel.save()
# artificial but better to test
# create a link table referencing rows in different spaces
ulabel.projects.add(project)


print("Created models")


# create a db connection url that works with RLS
def create_jwt_user(dsn_admin: str):
    db_role_handler = DbRoleHandler(dsn_admin)
    jwt_role_name = "permissions_jwt"
    jwt_db_url = db_role_handler.create(
        jwt_role_name, expires_in=None, alter_if_exists=True
    )
    db_role_handler.permission.grant_write_jwt(jwt_role_name)
    return jwt_db_url


pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url
jwt_db_url = create_jwt_user(pgurl)
ln.setup.settings.instance._db = jwt_db_url
ln.setup.settings.instance._persist()

print("Created jwt db connection")
