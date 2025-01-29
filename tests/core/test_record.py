import shutil
from inspect import signature
from pathlib import Path

import bionty as bt
import lamindb as ln
import pytest
from lamindb import _record
from lamindb._record import _search, suggest_records_with_similar_names
from lamindb.base.validation import FieldValidationError


def test_signatures():
    # this seems currently the easiest and most transparent
    # way to test violations of the signature equality
    # the MockORM class is needed to get inspect.signature
    # to work
    class Mock:
        pass

    # class methods
    class_methods = ["filter", "get", "df", "search", "lookup", "using"]
    for name in class_methods:
        setattr(Mock, name, getattr(_record, name))
        assert signature(getattr(Mock, name)) == _record.SIGS.pop(name)
    # methods
    for name, sig in _record.SIGS.items():
        assert signature(getattr(_record, name)) == sig


def test_validate_literal_fields():
    # validate literal
    with pytest.raises(FieldValidationError):
        ln.Transform(key="new-name-not-existing-123", type="invalid")


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


def test_search_and_get(get_search_test_filepaths):
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

    f = ln.Artifact.get(key="test-search5")
    f.suffix = ".txt"
    f.save()
    # multi-field search
    res = ln.Artifact.search("txt", field=["key", "description", "suffix"]).df()
    assert res.iloc[0].suffix == ".txt"

    # get

    artifact = ln.Artifact.get(description="test-search4")
    assert artifact == artifact4

    # because we're rendering Artifact.DoesNotExist private
    # in some use cases, we're not testing for it
    with pytest.raises(ln.Artifact._DoesNotExist):
        ln.Artifact.get(description="test-search1000000")

    #
    artifact0.delete(permanent=True, storage=True)
    artifact1.delete(permanent=True, storage=True)
    artifact2.delete(permanent=True, storage=True)
    artifact3.delete(permanent=True, storage=True)
    artifact4.delete(permanent=True, storage=True)
    artifact5.delete(permanent=True, storage=True)


def test_suggest_similar_names():
    ulabel1 = ln.ULabel(name="Test experiment 1").save()
    ulabel2 = ln.ULabel(name="Test experiment 2").save()
    ulabel3 = ln.ULabel(name="Special test experiment abc").save()
    ulabel4 = ln.ULabel(name="A very special test experiment abc").save()

    assert ln.ULabel(name="Test experiment 1").uid == ulabel1.uid

    assert suggest_records_with_similar_names(
        ulabel1, "name", {"name": "Test experiment 1"}
    )
    assert not suggest_records_with_similar_names(
        ulabel2, "name", {"name": "Test experiment 123"}
    )

    queryset = _search(
        ln.ULabel,
        "Test experiment 123",
        field="name",
        truncate_string=True,
        limit=3,
    )
    assert queryset.count() == 3

    queryset = _search(
        ln.ULabel,
        "Special test experiment abc",
        field="name",
        truncate_string=True,
        limit=3,
    )
    assert queryset.count() == 2
    assert queryset[0].name == "Special test experiment abc"

    ulabel1.delete()
    ulabel2.delete()
    ulabel3.delete()
    ulabel4.delete()


def test_pass_version():
    # creating a new transform on key retrieves the same transform
    # for as long as no source_code was saved
    transform = ln.Transform(key="mytransform", version="1").save()
    assert ln.Transform(key="mytransform", version="1") == transform
    # in case source code is saved
    transform.source_code = "dummy"
    transform.save()
    with pytest.raises(ValueError, match="Please increment the previous version"):
        ln.Transform(key="mytransform", version="1")


def test_get_name_field():
    transform = ln.Transform(key="test")
    transform.save()
    assert _record.get_name_field(ln.Run(transform)) == "started_at"
    with pytest.raises(ValueError):
        _record.get_name_field(ln.Artifact.ulabels.through())
    transform.delete()


def test_using():
    # the two below calls error if the records aren't found
    ln.Artifact.using("laminlabs/lamin-site-assets").get(1)
    ln.Artifact.using("laminlabs/lamin-site-assets").get(uid="MqEaGU7fXvxNy61R0000")
    # cross-database query
    hemangioblast = bt.CellType.from_source(name="hemangioblast").save()
    artifact = (
        ln.Artifact.using("laminlabs/lamin-dev")
        .filter(cell_types=hemangioblast)
        .first()
    )
    assert artifact is not None
    hemangioblast_dev = artifact.cell_types.get(name="hemangioblast")
    assert hemangioblast_dev.uid == hemangioblast.uid
    assert hemangioblast_dev.id != hemangioblast.id
    # query via list
    artifact_ref = (
        ln.Artifact.using("laminlabs/lamin-dev")
        .filter(cell_types__in=[hemangioblast])
        .first()
    )
    assert artifact == artifact_ref
