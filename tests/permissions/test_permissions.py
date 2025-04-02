import hubmodule.models as hm
import lamindb as ln
import psycopg2
import pytest
from django.db.utils import ProgrammingError
from jwt_utils import sign_jwt
from lamindb_setup.core.django import set_db_token

pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url
user_uuid = ln.setup.settings.user._uuid.hex
token = sign_jwt(pgurl, {"account_id": user_uuid})
set_db_token(token)


def test_fine_grained_permissions():
    # check select
    assert ln.ULabel.filter().count() == 3
    assert ln.Project.filter().count() == 2

    ulabel = ln.ULabel.get(name="space_all_ulabel")
    assert ulabel.projects.all().count() == 2
    # check delete
    # should delete
    ln.ULabel.get(name="full_access_ulabel").delete()
    assert ln.ULabel.filter().count() == 2
    # should not delete, does not error for some reason
    ln.ULabel.get(name="select_ulabel").delete()
    assert ln.ULabel.filter().count() == 2
    # space all
    ulabel.delete()
    assert ln.ULabel.filter().count() == 2
    # check insert
    # should succeed
    space = ln.models.Space.get(name="full access")
    ulabel = ln.ULabel(name="new label")
    ulabel.space = space
    ulabel.save()
    # should fail
    with pytest.raises(ProgrammingError):
        ln.ULabel(name="new label fail").save()
    for space_name in ["select access", "no access"]:
        space = ln.models.Space.get(name=space_name)
        ulabel = ln.ULabel(name="new label fail")
        ulabel.space = space
        with pytest.raises(ProgrammingError):
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
    with pytest.raises(ProgrammingError):
        ulabel.save()
    # space all
    ulabel = ln.ULabel.get(name="space_all_ulabel")
    ulabel.name = "space_all_ulabel update"
    with pytest.raises(ProgrammingError):
        ulabel.save()
    # check link tables
    # check insert
    project = ln.Project(name="Myproject")
    project.space = ln.models.Space.get(name="full access")
    project.save()
    ulabel = ln.ULabel.get(name="new label update")
    ulabel.projects.add(project)
    assert ulabel.projects.all().count() == 1
    # check select of a link table referencing unavailable rows
    assert ln.ULabel.get(name="select_ulabel").projects.all().count() == 0


def test_utility_tables():
    # can select in these tables
    assert ln.models.User.filter().count() == 1
    assert ln.models.Space.filter().count() == 4
    # can't select
    assert hm.Account.filter().count() == 0
    assert hm.Team.filter().count() == 0
    assert hm.AccountTeam.filter().count() == 0
    assert hm.AccessSpace.filter().count() == 0
    # can't update
    space = ln.models.Space.get(name="All")
    space.name = "new name"
    with pytest.raises(ProgrammingError):
        space.save()

    user = ln.models.User.filter().one()
    user.name = "new name"
    with pytest.raises(ProgrammingError):
        space.save()
    # can't delete
    ln.models.Space.get(name="All").delete()
    assert ln.models.Space.filter().count() == 4

    ln.models.User.filter().one().delete()
    assert ln.models.User.filter().count() == 1


def test_write_role():
    # switch user role to write
    with psycopg2.connect(pgurl) as conn, conn.cursor() as cur:
        cur.execute(
            "UPDATE hubmodule_account SET role = %s WHERE id = %s", ("write", user_uuid)
        )

    ln.ULabel(name="new label space all").save()
