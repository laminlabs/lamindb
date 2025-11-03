import lamindb as ln  # noqa
import hubmodule
import hubmodule.models as hm
from uuid import uuid4
from hubmodule._migrate import _apply_migrations_with_tracking
from hubmodule._setup import _setup_extensions, _setup_secret, _setup_utils_jwt
from hubmodule._rls import RLSGenerator
from laminhub_rest.core.postgres import DbRoleHandler
from pathlib import Path


# create a db connection url that works with RLS
instance_id = ln.setup.settings.instance._id


def create_jwt_user(dsn_admin: str, jwt_role_name: str):
    db_role_handler = DbRoleHandler(dsn_admin)
    jwt_db_url = db_role_handler.create(
        jwt_role_name, expires_in=None, alter_if_exists=True
    )
    db_role_handler.permission.grant_write_jwt(jwt_role_name)
    return jwt_db_url


pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url
jwt_role_name = f"{instance_id.hex}_jwt"
jwt_db_url = create_jwt_user(pgurl, jwt_role_name=jwt_role_name)

_setup_extensions(pgurl)
_setup_secret(pgurl)
_setup_utils_jwt(pgurl)
migrations_sql_dir = Path(hubmodule.__file__).parent / "sql/0004_migrations"
_apply_migrations_with_tracking(pgurl, migrations_sql_dir)

rls_generator = RLSGenerator(pgurl, jwt_role_name=jwt_role_name, public_role_name=None)

for i, table in enumerate(rls_generator._list_tables()):
    print(i, table.table_name, table.foreign_keys, table.has_space_id)

rls_generator.setup()


print("Created jwt db connection")

# create models

full_access = ln.Space(name="full access").save()  # type: ignore
select_access = ln.Space(name="select access").save()  # type: ignore
no_access = ln.Space(name="no access").save()  # type: ignore
# set read role for the default space
account = hm.Account(
    id=ln.setup.settings.user._uuid.hex, uid="accntid1", role="read"
).save()

# no access space
ulabel = ln.ULabel(name="no_access_ulabel")
ulabel.space = no_access
ulabel.save()
# set up access to this individual record with a dummy role,
# will work only after the role is changed to read, write or admin
hm.AccessRecord(
    account=account, record_type="lamindb_ulabel", record_id=ulabel.id, role="dummy"
).save()

project = ln.Project(name="No_access_project")  # type: ignore
project.space = no_access
project.save()

hm.AccessRecord(
    account=account, record_type="lamindb_project", record_id=project.id, role="dummy"
).save()

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

# default space, only select access by default
ulabel = ln.ULabel(name="default_space_ulabel").save()
ulabel.projects.add(project)

project = ln.Project(name="default_space_project").save()
ulabel.projects.add(project)

# create a link table referencing ulabel from the default space and project from select space
project = ln.Project(name="select_project")
project.space = select_access
project.save()

ulabel.projects.add(project)

# setup team and relevent models
team_access = ln.Space(name="team access").save()  # type: ignore
team = hm.Team(id=uuid4().hex, uid="teamuiduid11", name="test_team", role="read").save()
hm.AccountTeam(account=account, team=team).save()
hm.AccessSpace(team=team, space=team_access, role="read").save()

feature = ln.Feature(name="team_access_feature", dtype=float)
feature.space = team_access
feature.save()

# artifact for testing tracking error and artifactblock
artifact = ln.Artifact("README.md", description="test tracking error")
artifact.space = select_access
artifact.save()

# artifact for testing tracking error and locking
artifact = ln.Artifact(".gitignore", description="test locking")
artifact.space = full_access
artifact.is_locked = True
artifact.save()

# create a single record in the default space
record = ln.Record(name="test-record", is_type=False).save()
assert record.space_id == 1

print("Created models")

# save jwt db connection

ln.setup.settings.instance._db = jwt_db_url
ln.setup.settings.instance._persist()
