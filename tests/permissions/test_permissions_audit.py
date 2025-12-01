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
from lamindb.models.artifact import track_run_input
from lamindb_setup.core.django import DBToken, db_token_manager
from psycopg2.extensions import adapt

pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url

user_uuid = ln.setup.settings.user._uuid.hex
expiration = time.time() + 2000
# full collaborator token
token = sign_jwt(
    pgurl, {"account_id": user_uuid, "exp": expiration, "type": "collaborator"}
)
# read-only token
token_read = sign_jwt(
    pgurl, {"account_id": user_uuid, "exp": expiration, "type": "read-only"}
)
# init an instance of DBToken manually
db_token = DBToken({})
db_token._token = token
db_token._token_query = f"SELECT set_token({adapt(token).getquoted().decode()}, true);"
db_token._expiration = expiration

db_token_manager.set(db_token)


def test_authentication():
    # just check that the token was setup
    with connection.cursor() as cur:
        cur.execute(
            "SELECT 1 in (SELECT id FROM public.check_access() WHERE role = 'read');"
        )
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
    # test read-only token
    with connection.connection.cursor() as cur:
        cur.execute("SELECT set_token(%s); SELECT * FROM check_access()", (token_read,))
        result = cur.fetchall()
    assert len(result) == 1
    assert result[0] == (1, "read", "space")


def test_select_without_db_token():
    # with db token can be read in the default space
    with connection.cursor() as cur:
        cur.execute("SELECT * FROM lamindb_record;")
        results = cur.fetchall()
    assert len(results) == 1
    # the same
    assert ln.Record.filter().count() == 1
    # errors if can't select
    ln.Record.get(1)
    # no db token, everything in the default space
    with (
        pytest.raises(psycopg2.errors.RaiseException),
        connection.connection.cursor() as cur,
    ):
        cur.execute("SELECT * FROM lamindb_record;")
    with (
        pytest.raises(psycopg2.errors.RaiseException),
        connection.connection.cursor() as cur,
    ):
        cur.execute("SELECT * FROM lamindb_record WHERE id = 1;")
    # no db token, in different spaces
    with (
        pytest.raises(psycopg2.errors.RaiseException),
        connection.connection.cursor() as cur,
    ):
        cur.execute("SELECT * FROM lamindb_artifact;")
    with (
        pytest.raises(psycopg2.errors.RaiseException),
        connection.connection.cursor() as cur,
    ):
        cur.execute("SELECT * FROM lamindb_ulabel;")
    # no db token, utility tables
    with (
        pytest.raises(psycopg2.errors.RaiseException),
        connection.connection.cursor() as cur,
    ):
        cur.execute("SELECT * FROM lamindb_user;")
    with (
        pytest.raises(psycopg2.errors.RaiseException),
        connection.connection.cursor() as cur,
    ):
        cur.execute("SELECT * FROM lamindb_space;")


def test_fine_grained_permissions_account():
    # check select
    assert ln.ULabel.filter().count() == 3
    assert ln.Project.filter().count() == 2

    ulabel = ln.ULabel.get(name="default_space_ulabel")
    assert ulabel.projects.all().count() == 2
    # check delete
    # should delete
    ln.ULabel.get(name="full_access_ulabel").delete(permanent=True)
    assert ln.ULabel.filter().count() == 2
    # should not delete, does not error for some reason
    ln.ULabel.get(name="select_ulabel").delete(permanent=True)
    assert ln.ULabel.filter().count() == 2
    # default space
    ulabel.delete(permanent=True)
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
    # test RootBlock, can do due to write access to some spaces
    root_block = ln.models.RootBlock(name="instance", content="test").save()
    root_block.content = "test 2"
    root_block.save()
    # test SpaceBlock
    space = ln.Space.get(name="select access")
    with pytest.raises(ln.errors.NoWriteAccess):
        ln.models.SpaceBlock(space=space, content="test").save()
    # test ArtifactBlock, artifact is read-only
    artifact = ln.Artifact.get(description="test tracking error")
    with pytest.raises(ProgrammingError):
        ln.models.ArtifactBlock(artifact=artifact, content="test").save()
    # test BranchBlock, the account is read-only
    branch = ln.Branch.get(1)  # main branch in all space
    with pytest.raises(ProgrammingError):
        ln.models.BranchBlock(branch=branch, content="test").save()


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

    ulabel.delete(permanent=True)
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
        ln.Space(name="new space").save()

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

    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            "UPDATE hubmodule_team SET role = 'read' WHERE uid = 'teamuiduid11'",
        )


def test_locking():
    artifact = ln.Artifact.get(description="test locking")
    artifact.description = "new description"
    with pytest.raises(ln.errors.NoWriteAccess) as e:
        artifact.save()
    assert "It is not allowed to modify or create locked" in str(e)


def test_tracking_error():
    # switch user role to write to create the transform and run
    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            "UPDATE hubmodule_account SET role = 'write' WHERE id = %s", (user_uuid,)
        )

    artifact = ln.Artifact.get(description="test tracking error")

    transform = ln.Transform(key="My transform").save()
    run = ln.Run(transform).save()

    # this error because ln.setup.settings.instance._db_permissions is not jwt
    # it is None
    with pytest.raises(ln.errors.NoWriteAccess) as e:
        track_run_input(artifact, run)
    assert "You’re not allowed to write to the instance " in str(e)

    # the instance is local so we set this manually
    ln.setup.settings.instance._db_permissions = "jwt"
    # artifact.space is not available for writes
    with pytest.raises(ln.errors.NoWriteAccess) as e:
        track_run_input(artifact, run)
    assert "You’re not allowed to write to the space " in str(e)

    # this artifact is locked
    artifact = ln.Artifact.get(description="test locking")
    with pytest.raises(ln.errors.NoWriteAccess) as e:
        track_run_input(artifact, run)
    assert "It is not allowed to modify locked records" in str(e)

    # switch user role back to read
    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            "UPDATE hubmodule_account SET role = 'read' WHERE id = %s", (user_uuid,)
        )
    # as the user is read-only now, 2 spaces are unavailable for writes (artifact.space, run.space)
    artifact = ln.Artifact.get(description="test tracking error")
    with pytest.raises(ln.errors.NoWriteAccess) as e:
        track_run_input(artifact, run)
    assert "You’re not allowed to write to the spaces " in str(e)

    ln.setup.settings.instance._db_permissions = None


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
    subprocess.run(  # noqa: S602
        f"python {script1_path}",
        shell=True,
        check=True,
    )
