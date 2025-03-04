import lamindb as ln
import pytest
from django.db.utils import ProgrammingError
from jwt_utils import sign_jwt
from lamindb_setup.core.django import set_token

pgurl = "postgresql://postgres:pwd@0.0.0.0:5432/pgtest"  # admin db connection url
token = sign_jwt(pgurl, {"account_id": ln.setup.settings.user._uuid.hex})
set_token(token)


def test_fine_grained_permissions():
    # check select
    assert ln.ULabel.filter().count() == 2
    assert ln.Project.filter().count() == 0
    # check delete
    # should delete
    ln.ULabel.get(name="full_access_ulabel").delete()
    assert ln.ULabel.filter().count() == 1
    # should not delete, does not error for some reason
    ln.ULabel.get(name="select_ulabel").delete()
    assert ln.ULabel.filter().count() == 1
    # check insert
    # should succeed
    space = ln.models.Space.get(name="full access")
    ulabel = ln.ULabel(name="new label")
    ulabel.space = space
    ulabel.save()
    # should fail
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
