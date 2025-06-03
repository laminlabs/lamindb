import subprocess
import time
from pathlib import Path
from uuid import uuid4

import hubmodule.models as hm
import lamindb as ln
import psycopg2
import pytest
from django.db import connection, transaction
from django.db.utils import InternalError, ProgrammingError
from jwt_utils import sign_jwt
from lamindb_setup.core.django import DBToken, db_token_manager
from psycopg2.extensions import adapt

pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url

user_uuid = ln.setup.settings.user._uuid.hex
expiration = time.time() + 2000
token = sign_jwt(pgurl, {"account_id": user_uuid, "exp": expiration})
# init an instance of DBToken manually
db_token = DBToken({})
db_token._token = token
db_token._token_query = f"SELECT set_token({adapt(token).getquoted().decode()}, true);"
db_token._expiration = expiration

db_token_manager.set(db_token)


def test_authentication():
    # just check that the token was setup
    with connection.cursor() as cur:
        cur.execute("SELECT get_account_id();")
        account_id = cur.fetchall()[0][0]
    assert account_id.hex == user_uuid
    # test that auth can't be hijacked
    # false table created before
    with (
        pytest.raises(psycopg2.errors.DuplicateTable),
        connection.connection.cursor() as cur,
    ):
        cur.execute(
            """
            CREATE TEMP TABLE account_id(val uuid PRIMARY KEY) ON COMMIT DROP;
            SELECT set_token(%s);
            """,
            (token,),
        )
    # check that jwt user can't set arbitrary account_id manually
    with (
        pytest.raises(psycopg2.errors.RaiseException),
        connection.connection.cursor() as cur,
    ):
        cur.execute(
            """
            CREATE TEMP TABLE account_id(val uuid PRIMARY KEY) ON COMMIT DROP;
            INSERT INTO account_id(val) VALUES (gen_random_uuid());
            SELECT get_account_id();
            """
        )
    # check manual insert
    with (
        pytest.raises(psycopg2.errors.InsufficientPrivilege),
        connection.connection.cursor() as cur,
    ):
        cur.execute(
            """
            SELECT set_token(%s);
            INSERT INTO account_id(val) VALUES (gen_random_uuid());
            """,
            (token,),
        )
    # test access to the security schema
    with (
        pytest.raises(psycopg2.errors.InsufficientPrivilege),
        connection.connection.cursor() as cur,
    ):
        cur.execute("SELECT security.get_secret('jwt_secret');")


def test_fine_grained_permissions_account():
    # check select
    assert ln.ULabel.filter().count() == 3
    assert ln.Project.filter().count() == 2

    ulabel = ln.ULabel.get(name="default_space_ulabel")
    assert ulabel.projects.all().count() == 2
    # check delete
    # should delete
    ln.ULabel.get(name="full_access_ulabel").delete()
    assert ln.ULabel.filter().count() == 2
    # should not delete, does not error for some reason
    ln.ULabel.get(name="select_ulabel").delete()
    assert ln.ULabel.filter().count() == 2
    # default space
    ulabel.delete()
    assert ln.ULabel.filter().count() == 2
    # check insert
    # should succeed
    space = ln.Space.get(name="full access")
    ulabel = ln.ULabel(name="new label")
    ulabel.space = space
    ulabel.save()
    # should fail
    with pytest.raises(ln.errors.NoWriteAccess):
        ln.ULabel(name="new label fail").save()
    for space_name in ["select access", "no access"]:
        space = ln.Space.get(name=space_name)
        ulabel = ln.ULabel(name="new label fail")
        ulabel.space = space
        with pytest.raises(ln.errors.NoWriteAccess):
            ulabel.save()
    # check update
    # should succeed
    ulabel = ln.ULabel.get(name="new label")
    ulabel.name = "new label update"
    ulabel.save()
    ulabel = ln.ULabel.get(name="new label update")  # check that it is saved
    # should fail
    ulabel = ln.ULabel.get(name="select_ulabel")
    ulabel.name = "select_ulabel update"
    with pytest.raises(ln.errors.NoWriteAccess):
        ulabel.save()
    # default space
    ulabel = ln.ULabel.get(name="default_space_ulabel")
    ulabel.name = "default_space_ulabel update"
    with pytest.raises(ln.errors.NoWriteAccess):
        ulabel.save()
    # check link tables
    # check insert
    project = ln.Project(name="Myproject")
    project.space = ln.Space.get(name="full access")
    project.save()
    ulabel = ln.ULabel.get(name="new label update")
    ulabel.projects.add(project)
    assert ulabel.projects.all().count() == 1
    # check select of a link table referencing unavailable rows
    assert ln.ULabel.get(name="select_ulabel").projects.all().count() == 0


def test_fine_grained_permissions_team():
    assert ln.Feature.filter().count() == 1
    ln.Feature.get(name="team_access_feature")


# tests that token is set properly in atomic blocks
def test_atomic():
    with transaction.atomic():
        assert ln.Feature.filter().count() == 1
        # test with nested
        with transaction.atomic():
            assert ln.Feature.filter().count() == 1

            feature = ln.Feature(name="atomic_feature", dtype=float)
            feature.space = ln.Space.get(name="full access")
            feature.save()

    assert ln.Feature.filter().count() == 2


def test_utility_tables():
    # can select in these tables
    assert ln.models.User.filter().count() == 1
    assert ln.Space.filter().count() == 5
    # can't select
    assert hm.Account.filter().count() == 0
    assert hm.Team.filter().count() == 0
    assert hm.AccountTeam.filter().count() == 0
    assert hm.AccessSpace.filter().count() == 0
    # can't update
    space = ln.Space.get(id=1)  # default space
    space.name = "new name"
    with pytest.raises(ProgrammingError):
        space.save()

    user = ln.models.User.filter().one()
    user.name = "new name"
    with pytest.raises(ProgrammingError):
        space.save()
    # can't insert
    with pytest.raises(ProgrammingError):
        ln.Space(name="new space", uid="00000005").save()

    with pytest.raises(ProgrammingError):
        hm.Account(id=uuid4().hex, uid="accntid2", role="admin").save()


def test_write_role():
    # switch user role to write
    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            "UPDATE hubmodule_account SET role = 'write' WHERE id = %s", (user_uuid,)
        )

    ln.ULabel(name="new label account default space").save()

    # switch user role back to read and team role to write
    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            "UPDATE hubmodule_account SET role = 'read' WHERE id = %s", (user_uuid,)
        )
        cur.execute(
            "UPDATE hubmodule_team SET role = 'write' WHERE uid = 'teamuiduid11'",
        )

    ln.ULabel(name="new label team default space").save()


def test_token_reset():
    db_token_manager.reset()

    # account_id is not set
    with pytest.raises(InternalError) as error:
        ln.ULabel.filter().count()
    assert "JWT is not set" in error.exconly()

    with pytest.raises(InternalError) as error, transaction.atomic():
        ln.ULabel.filter().count()
    assert "JWT is not set" in error.exconly()


def test_lamin_dev():
    script1_path = Path(__file__).parent.resolve() / "scripts/check_lamin_dev.py"
    script2_path = Path(__file__).parent.resolve() / "scripts/clean_lamin_dev.py"
    # TODO: if we don't access the instance here, it will be changed
    subprocess.run(  # noqa: S602
        f"python {script1_path}",
        shell=True,
        check=True,
    )
    result = subprocess.run(  # noqa: S602
        "lamin save .gitignore --key mytest --space 'Our test space for CI'",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert "key='mytest'" in result.stdout.decode()
    assert "storage path:" in result.stdout.decode()
    assert result.returncode == 0

    result = subprocess.run(  # noqa: S602
        f"python {script2_path}",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 0
