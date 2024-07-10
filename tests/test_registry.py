import shutil
from inspect import signature
from pathlib import Path

import lamindb as ln
import pytest
from lamindb import _registry as registry


def test_signatures():
    # this seems currently the easiest and most transparent
    # way to test violations of the signature equality
    # the MockORM class is needed to get inspect.signature
    # to work
    class Mock:
        pass

    # class methods
    class_methods = ["search", "lookup", "from_values", "using"]
    for name in class_methods:
        setattr(Mock, name, getattr(registry, name))
        assert signature(getattr(Mock, name)) == registry.SIGS.pop(name)
    # methods
    for name, sig in registry.SIGS.items():
        assert signature(getattr(registry, name)) == sig


def test_init_with_args():
    with pytest.raises(ValueError) as error:
        # can't use ULabel here because it raises "Only one non-keyword arg allowed"
        ln.User("an arg")
    assert (
        error.exconly()
        == "ValueError: please provide keyword arguments, not plain arguments"
    )


def test_validate_required_fields():
    # label has a required name
    with pytest.raises(TypeError):
        ln.ULabel()
    # label has a required name
    with pytest.raises(TypeError):
        ln.ULabel(description="test")


@pytest.fixture
def get_search_test_filepaths():
    Path("unregistered_storage/").mkdir(exist_ok=True)
    filepaths = [Path(f"./unregistered_storage/test-search{i}") for i in range(6)]
    for filepath in filepaths:
        filepath.write_text(filepath.name)
    yield None
    shutil.rmtree("unregistered_storage/")


def test_search_artifact(get_search_test_filepaths):
    artifact1 = ln.Artifact(
        "./unregistered_storage/test-search1", description="nonsense"
    )
    artifact1.save()
    artifact2 = ln.Artifact(
        "./unregistered_storage/test-search2", description="nonsense"
    )
    artifact2.save()

    # on purpose to be search3 to test duplicated search
    artifact0 = ln.Artifact(
        "./unregistered_storage/test-search0", description="test-search3"
    )
    artifact0.save()
    artifact3 = ln.Artifact(
        "./unregistered_storage/test-search3", description="test-search3"
    )
    artifact3.save()
    artifact4 = ln.Artifact(
        "./unregistered_storage/test-search4", description="test-search4"
    )
    artifact4.save()

    result = ln.Artifact.search("search3").df()
    assert result.iloc[0].description == "test-search3"
    assert result.iloc[1].description == "test-search3"

    # no returning entries if all search results have __ratio__ 0
    # need a better search string below
    # assert ln.Artifact.search("x").shape[0] == 0

    artifact5 = ln.Artifact("./unregistered_storage/test-search5", key="test-search5")
    artifact5.save()
    res = ln.Artifact.search("search5").df()
    assert res.iloc[0].key == "test-search5"

    res_q = ln.Artifact.search("search5")
    assert res_q[0].key == "test-search5"
    # queryset returns the same order of results
    assert res.uid.tolist() == [i.uid for i in res_q]

    f = ln.Artifact.filter(key="test-search5").one()
    f.suffix = ".txt"
    f.save()
    # multi-field search
    res = ln.Artifact.search("txt", field=["key", "description", "suffix"]).df()
    assert res.iloc[0].suffix == ".txt"
    artifact0.delete(permanent=True, storage=True)
    artifact1.delete(permanent=True, storage=True)
    artifact2.delete(permanent=True, storage=True)
    artifact3.delete(permanent=True, storage=True)
    artifact4.delete(permanent=True, storage=True)
    artifact5.delete(permanent=True, storage=True)


def test_pass_version():
    transform = ln.Transform(name="mytransform", version="1")
    transform.save()
    assert ln.Transform(name="mytransform", version="1") == transform


def test_get_name_field():
    transform = ln.Transform(name="test")
    transform.save()
    assert registry.get_name_field(ln.Run(transform)) == "created_at"
    with pytest.raises(ValueError):
        registry.get_name_field(ln.Artifact.ulabels.through())
    transform.delete()
