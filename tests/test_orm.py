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
        ln.Label("an arg", name="test")


def test_validate_required_fields():
    # label has a required name
    with pytest.raises(TypeError):
        ln.Label()
    # label has a required name
    with pytest.raises(TypeError):
        ln.Label(description="test")


def test_search_file():
    for i in range(6):
        with open(f"test-search{i}", "w") as f:
            f.write(f"test-search{i}")

    file1 = ln.File("test-search1")
    file1.save()
    file2 = ln.File("test-search2")
    file2.save()

    # on purpose to be search3 to test duplicated search
    file0 = ln.File("test-search0", description="test-search3")
    file0.save()
    file3 = ln.File("test-search3", description="test-search3")
    file3.save()
    file4 = ln.File("test-search4", description="test-search4")
    file4.save()

    res = ln.File.search("search3")
    assert res.iloc[0].description == "test-search3"
    assert res.iloc[1].description == "test-search3"

    # no returning entries if all search results have __ratio__ 0
    assert ln.File.search("x").shape[0] == 0

    file5 = ln.File("test-search5", key="test-search5")
    file5.save()
    res = ln.File.search("search5")
    assert res.iloc[0].key == "test-search5"

    res_q = ln.File.search("search5", return_queryset=True)
    assert res_q[0].key == "test-search5"
    # queryset returns the same order of results
    assert res.index.tolist() == [i.id for i in res_q]

    f = ln.File.select(key="test-search5").one()
    f.suffix = ".txt"
    f.save()
    # multi-field search
    res = ln.File.search("txt", field=["key", "description", "suffix"])
    assert res.iloc[0].suffix == ".txt"
    file0.delete(storage=True)
    file1.delete(storage=True)
    file2.delete(storage=True)
    file3.delete(storage=True)
    file4.delete(storage=True)
    file5.delete(storage=True)
