from inspect import signature

import pytest

import lamindb as ln
from lamindb import _orm as orm


def test_signatures():
    # this seems currently the easiest and most transparent
    # way to test violations of the signature equality
    # the MockORM class is needed to get inspect.signature
    # to work
    class Mock:
        pass

    # class methods
    class_methods = ["search", "lookup", "from_values", "inspect", "map_synonyms"]
    for name in class_methods:
        setattr(Mock, name, getattr(orm, name))
        assert signature(getattr(Mock, name)) == orm.SIGS.pop(name)
    # methods
    for name, sig in orm.SIGS.items():
        assert signature(getattr(orm, name)) == sig


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
