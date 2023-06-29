from inspect import signature

import pytest

import lamindb as ln
from lamindb import _orm as orm


def test_signatures():
    # this seems currently the easiest and most transparent
    # way to test violations of the signature equality
    # the MockORM class is needed to get inspect.signature
    # to work
    class MockORM:
        pass

    MockORM.search = orm.search
    assert signature(MockORM.search) == orm.SIG_ORM_SEARCH
    MockORM.lookup = orm.lookup
    assert signature(MockORM.lookup) == orm.SIG_ORM_LOOKUP
    MockORM.from_values = orm.from_values
    assert signature(MockORM.from_values) == orm.SIG_ORM_FROM_VALUES
    MockORM.add_synonym = orm.from_values
    assert signature(MockORM.add_synonym) == orm.SIG_ORM_ADD_SYNONYM
    MockORM.remove_synonym = orm.remove_synonym
    assert signature(MockORM.remove_synonym) == orm.SIG_ORM_REMOVE_SYNONYM


def test_init_with_args():
    with pytest.raises(ValueError):
        ln.Tag("an arg")


def test_validate_required_fields():
    # tag has a required name
    with pytest.raises(TypeError):
        ln.Tag()
    # project has a required name
    with pytest.raises(TypeError):
        ln.Project(external_id="test")
