import re
import shutil
from pathlib import Path

import bionty as bt
import lamindb as ln
import pytest
from lamindb.errors import FieldValidationError
from lamindb.models.sqlrecord import (
    _get_record_kwargs,
    _search,
    get_name_field,
    suggest_records_with_similar_names,
)


def test_validate_literal_fields():
    # validate literal
    with pytest.raises(FieldValidationError):
        ln.Transform(key="new-name-not-existing-123", type="invalid")


def test_init_with_args():
    with pytest.raises(
        FieldValidationError,
        match=re.escape(
            "Use keyword arguments instead of positional arguments, e.g.: User(name='...')"
        )
        + r".*",
    ):
        # can't use ULabel here because it raises "Only one non-keyword arg allowed"
        ln.User("an arg")


def test_validate_required_fields():
    # label has a required name
    with pytest.raises(FieldValidationError):
        ln.ULabel()
    # label has a required name
    with pytest.raises(FieldValidationError):
        ln.ULabel(description="test")


@pytest.fixture
def get_search_test_filepaths():
    Path("unregistered_storage/").mkdir(exist_ok=True)
    filepaths = [Path(f"./unregistered_storage/test-search{i}.txt") for i in range(6)]
    for filepath in filepaths:
        filepath.write_text(filepath.name)
    yield None
    shutil.rmtree("unregistered_storage/")


def test_search_and_get(get_search_test_filepaths):
    artifact1 = ln.Artifact(
        "./unregistered_storage/test-search1.txt", description="nonsense"
    )
    artifact1.save()
    artifact2 = ln.Artifact(
        "./unregistered_storage/test-search2.txt", description="nonsense"
    )
    artifact2.save()

    # on purpose to be search3 to test duplicated search
    artifact0 = ln.Artifact(
        "./unregistered_storage/test-search0.txt", description="test-search3"
    )
    artifact0.save()
    artifact3 = ln.Artifact(
        "./unregistered_storage/test-search3.txt", description="test-search3"
    )
    artifact3.save()
    artifact4 = ln.Artifact(
        "./unregistered_storage/test-search4.txt", description="test-search4"
    )
    artifact4.save()

    result = ln.Artifact.search("search3").df()
    assert result.iloc[0].description == "test-search3"
    assert result.iloc[1].description == "test-search3"

    # no returning entries if all search results have __ratio__ 0
    # need a better search string below
    # assert ln.Artifact.search("x").shape[0] == 0

    artifact5 = ln.Artifact(
        "./unregistered_storage/test-search5.txt", key="test-search5.txt"
    )
    artifact5.save()
    res = ln.Artifact.search("search5").df()
    assert res.iloc[0].key == "test-search5.txt"

    res_q = ln.Artifact.search("search5")
    assert res_q[0].key == "test-search5.txt"
    # queryset returns the same order of results
    assert res.uid.tolist() == [i.uid for i in res_q]

    # multi-field search
    res = ln.Artifact.search("txt", field=["key", "description", "suffix"]).df()
    assert res.iloc[0].suffix == ".txt"

    # get

    artifact = ln.Artifact.get(description="test-search4")
    assert artifact == artifact4

    with pytest.raises(ln.Artifact.DoesNotExist):
        ln.Artifact.get(description="test-does-not-exist")

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
    transform = ln.Transform(key="test").save()
    assert get_name_field(ln.Run(transform)) == "started_at"
    with pytest.raises(ValueError):
        get_name_field(ln.Artifact.ulabels.through())
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
    # check that .using provided with the current intance does nothing
    assert ln.User.using("lamindb-unit-tests-core").first()._state.db == "default"
    user = ln.setup.settings.user.handle
    assert (
        ln.User.using(f"{user}/lamindb-unit-tests-core").first()._state.db == "default"
    )


def test_get_record_kwargs():
    assert _get_record_kwargs(ln.Feature) == [
        ("name", "str"),
        ("dtype", "Dtype | Registry | list[Registry] | FieldAttr"),
        ("type", "Feature | None"),
        ("is_type", "bool"),
        ("unit", "str | None"),
        ("description", "str | None"),
        ("synonyms", "str | None"),
        ("nullable", "bool"),
        (
            "default_value",
            "str | None",
        ),
        ("coerce_dtype", "bool"),
        (
            "cat_filters",
            "dict[str",
        ),
    ]


def test_get_record_kwargs_empty():
    class EmptySQLRecord:
        pass

    assert _get_record_kwargs(EmptySQLRecord) == []

    class NoInitSQLRecord:
        def method(self):
            pass

    assert _get_record_kwargs(NoInitSQLRecord) == []
