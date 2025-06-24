import subprocess
import time
from pathlib import Path
from uuid import uuid4

import hubmodule.models as hm
import lamindb as ln
import psycopg2
import pytest
from django.db import connection, transaction
from django.db.utils import IntegrityError, InternalError, ProgrammingError
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
        cur.execute("SELECT 1 in (SELECT id FROM check_access() WHERE role = 'read');")
        result = cur.fetchall()[0][0]
    assert result
    # check querying without setting jwt
    with (
        pytest.raises(psycopg2.errors.RaiseException),
        connection.connection.cursor() as cur,
    ):
        cur.execute("SELECT * FROM lamindb_ulabel;")
    # test that auth can't be hijacked
    # false table created before
    with (
        pytest.raises(psycopg2.errors.DuplicateTable),
        connection.connection.cursor() as cur,
    ):
        cur.execute(
            """
            CREATE TEMP TABLE access(
                id int,
                role varchar(20),
                type text
            ) ON COMMIT DROP;
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
            CREATE TEMP TABLE access(
                id int,
                role varchar(20),
                type text
            ) ON COMMIT DROP;
            INSERT INTO access (id, role, type)
            VALUES (1, 'admin', 'space');
            SELECT * FROM check_access();
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
            INSERT INTO access (id, role, type)
            VALUES (1, 'admin', 'space');
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


def test_fine_grained_permissions_single_records():
    assert not ln.ULabel.filter(name="no_access_ulabel").exists()
    assert not ln.Project.filter(name="No_access_project").exists()

    # switch access to this ulabel to read
    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            """
            UPDATE hubmodule_accessrecord SET role = 'read'
            WHERE account_id = %s AND record_type = 'lamindb_ulabel'
            """,
            (user_uuid,),
        )

    ulabel = ln.ULabel.get(name="no_access_ulabel")

    new_name = "new_name_single_rls_access_ulabel"
    ulabel.name = new_name
    with pytest.raises(ln.errors.NoWriteAccess):
        ulabel.save()

    # switch access to this ulabel to write
    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            """
            UPDATE hubmodule_accessrecord SET role = 'write'
            WHERE account_id = %s AND record_type = 'lamindb_ulabel'
            """,
            (user_uuid,),
        )

    ulabel.save()

    # switch access to this ulabel to write
    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            """
            UPDATE hubmodule_accessrecord SET role = 'read'
            WHERE account_id = %s AND record_type = 'lamindb_project'
            """,
            (user_uuid,),
        )

    project = ln.Project.get(name="No_access_project")
    # can't insert into lamindb_ulabelproject because the project is still read-only
    with pytest.raises(ProgrammingError):
        ulabel.projects.add(project)

    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            """
            UPDATE hubmodule_accessrecord SET role = 'write'
            WHERE account_id = %s AND record_type = 'lamindb_project'
            """,
            (user_uuid,),
        )

    ulabel.projects.add(project)
    assert ulabel.projects.count() == 1

    ulabel.delete()
    assert not ln.ULabel.filter(name="no_access_ulabel").exists()


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
    assert ln.User.filter().count() == 1
    assert ln.Space.filter().count() == 5
    # can't select
    assert hm.Account.filter().count() == 0
    assert hm.Team.filter().count() == 0
    assert hm.AccountTeam.filter().count() == 0
    assert hm.AccessSpace.filter().count() == 0
    assert hm.AccessRecord.filter().count() == 0
    # can't update a space
    space = ln.Space.get(id=1)  # default space
    space.name = "new name"
    with pytest.raises(ProgrammingError):
        space.save()
    # can't update a user
    user = ln.User.filter().one()
    user.name = "new name"
    # as we allow insert but not update on the user table
    # it looks like the db raises IntegrityError insead of the rls error
    # because just tries to insert with the same id and fails
    with pytest.raises(IntegrityError):
        user.save()
    # can insert a user because has write access to a space
    ln.User(handle="insert_new_user", uid="someuidd").save()
    assert ln.User.filter().count() == 2
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
    # connect to the instance before saving
    subprocess.run(  # noqa: S602
        "lamin connect laminlabs/lamin-dev",
        shell=True,
        check=True,
    )
    result = subprocess.run(  # noqa: S602
        "lamin save .gitignore --key mytest --space 'Our test space for CI'",
        shell=True,
        capture_output=True,
    )
    assert "key='mytest'" in result.stdout.decode()
    assert "storage path:" in result.stdout.decode()
    assert result.returncode == 0

    result = subprocess.run(  # noqa: S602
        f"python {script2_path}",
        shell=True,
        capture_output=False,
    )
    assert result.returncode == 0
