import re
import shutil
import textwrap
from pathlib import Path

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.errors import FieldValidationError
from lamindb.models.sqlrecord import (
    _get_record_kwargs,
    _search,
    get_name_field,
    suggest_records_with_similar_names,
)


def test_feature_describe():
    description = textwrap.dedent("""\
    Feature
      Simple fields
        .uid: CharField
        .name: CharField
        .unit: CharField
        .description: TextField
        .array_rank: SmallIntegerField
        .array_size: IntegerField
        .array_shape: JSONField
        .synonyms: TextField
        .default_value: JSONField
        .nullable: BooleanField
        .coerce: BooleanField
        .is_type: BooleanField
        .is_locked: BooleanField
        .created_at: DateTimeField
        .updated_at: DateTimeField
      Relational fields
        .branch: Branch
        .space: Space
        .created_by: User
        .run: Run
        .type: Feature
        .schemas: Schema
        .features: Feature
        .values: JsonValue
        .projects: Project
        .ablocks: FeatureBlock
    """).strip()
    assert description == ln.Feature.describe(return_str=True)


def test_artifact_describe():
    description = textwrap.dedent("""\
    Artifact
      Simple fields
        .uid: CharField
        .key: CharField
        .description: TextField
        .suffix: CharField
        .kind: CharField
        .otype: CharField
        .size: BigIntegerField
        .hash: CharField
        .n_files: BigIntegerField
        .n_observations: BigIntegerField
        .version_tag: CharField
        .is_latest: BooleanField
        .is_locked: BooleanField
        .created_at: DateTimeField
        .updated_at: DateTimeField
      Relational fields
        .branch: Branch
        .space: Space
        .storage: Storage
        .run: Run
        .schema: Schema
        .created_by: User
        .input_of_runs: Run
        .recreating_runs: Run
        .schemas: Schema
        .json_values: JsonValue
        .artifacts: Artifact
        .linked_in_records: Record
        .users: User
        .runs: Run
        .ulabels: ULabel
        .linked_by_artifacts: Artifact
        .collections: Collection
        .records: Record
        .references: Reference
        .projects: Project
        .ablocks: ArtifactBlock
      Bionty fields
        .organisms: bionty.Organism
        .genes: bionty.Gene
        .proteins: bionty.Protein
        .cell_markers: bionty.CellMarker
        .tissues: bionty.Tissue
        .cell_types: bionty.CellType
        .diseases: bionty.Disease
        .cell_lines: bionty.CellLine
        .phenotypes: bionty.Phenotype
        .pathways: bionty.Pathway
        .experimental_factors: bionty.ExperimentalFactor
        .developmental_stages: bionty.DevelopmentalStage
        .ethnicities: bionty.Ethnicity
    """).strip()
    assert description == ln.Artifact.describe(return_str=True)


def test_repr_describe():
    user = ln.User.filter().first()
    assert user.__repr__().startswith("User")
    assert user.describe(return_str=True).startswith("User")


def test_validate_literal_fields():
    # validate literal
    with pytest.raises(FieldValidationError):
        ln.Transform(key="new-name-not-existing-123", kind="invalid")


def test_init_with_args():
    with pytest.raises(
        FieldValidationError,
        match=re.escape(
            "Use keyword arguments instead of positional arguments, e.g.: User(name='...')"
        )
        + r".*",
    ):
        # can't use Record here because it raises "Only one non-keyword arg allowed"
        ln.User("an arg")


def test_validate_required_fields():
    # ULabel has a required name
    with pytest.raises(FieldValidationError):
        ln.ULabel()
    # ULabel has a required name
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

    result = ln.Artifact.search("search3").to_dataframe()
    assert result.iloc[0].description == "test-search3"
    assert result.iloc[1].description == "test-search3"

    # no returning entries if all search results have __ratio__ 0
    # need a better search string below
    # assert ln.Artifact.search("x").shape[0] == 0

    artifact5 = ln.Artifact(
        "./unregistered_storage/test-search5.txt", key="test-search5.txt"
    )
    artifact5.save()
    res = ln.Artifact.search("search5").to_dataframe()
    assert res.iloc[0].key == "test-search5.txt"

    res_q = ln.Artifact.search("search5")
    assert res_q[0].key == "test-search5.txt"
    # queryset returns the same order of results
    assert res.uid.tolist() == [i.uid for i in res_q]

    # multi-field search
    res = ln.Artifact.search(
        "txt", field=["key", "description", "suffix"]
    ).to_dataframe()
    assert res.iloc[0].suffix == ".txt"

    # get

    artifact = ln.Artifact.get(description="test-search4")
    assert artifact == artifact4

    with pytest.raises(ln.errors.ObjectDoesNotExist):
        ln.Artifact.get(description="test-does-not-exist")

    artifact0.delete(permanent=True, storage=True)
    artifact1.delete(permanent=True, storage=True)
    artifact2.delete(permanent=True, storage=True)
    artifact3.delete(permanent=True, storage=True)
    artifact4.delete(permanent=True, storage=True)
    artifact5.delete(permanent=True, storage=True)


def test_suggest_similar_names():
    record1 = ln.Record(name="Test experiment 1").save()
    record2 = ln.Record(name="Test experiment 2").save()
    record3 = ln.Record(name="Special test experiment abc").save()
    record4 = ln.Record(name="A very special test experiment abc").save()

    assert ln.Record(name="Test experiment 1").uid == record1.uid

    assert suggest_records_with_similar_names(
        record1, "name", {"name": "Test experiment 1"}
    )
    assert not suggest_records_with_similar_names(
        record2, "name", {"name": "Test experiment 123"}
    )

    queryset = _search(
        ln.Record,
        "Test experiment 123",
        field="name",
        truncate_string=True,
        limit=3,
    )
    assert queryset.count() == 3

    queryset = _search(
        ln.Record,
        "Special test experiment abc",
        field="name",
        truncate_string=True,
        limit=3,
    )
    assert queryset.count() == 2
    assert queryset[0].name == "Special test experiment abc"

    record1.delete(permanent=True)
    record2.delete(permanent=True)
    record3.delete(permanent=True)
    record4.delete(permanent=True)


def test_pass_version():
    # creating a new transform on key retrieves the same transform
    # for as long as no source_code was saved
    transform = ln.Transform(key="mytransform", version="1").save()
    assert transform.version_tag == "1"
    assert transform.version == "1"
    assert ln.Transform(key="mytransform", version="1") == transform
    # in case source code is saved
    transform.source_code = "dummy"
    transform.save()
    with pytest.raises(ValueError) as e:
        ln.Transform(key="mytransform", version="1")
    assert (
        e.exconly()
        == "ValueError: Please change the version tag or leave it `None`, '1' is already taken"
    )


def test_delete():
    record = ln.Record(name="test-delete")
    # record not yet saved, delete has no effect
    result = record.delete()
    assert result is None
    assert record.branch_id == 1
    record.save()
    result = record.delete()
    assert result is None
    assert record.branch_id == -1
    result = record.delete(permanent=True)
    assert isinstance(result, tuple)
    assert len(result) == 2
    deleted_count, deleted_dict = result
    assert deleted_count == 1
    assert isinstance(deleted_dict, dict)
    assert ln.Record.filter(name="test-delete").exists() is False


def test_get_name_field():
    transform = ln.Transform(key="test").save()
    assert get_name_field(ln.Run(transform)) == "started_at"
    with pytest.raises(ValueError):
        get_name_field(ln.Artifact.records.through())
    transform.delete(permanent=True)


def test_using():
    # the two below calls error if the records aren't found
    ln.Artifact.connect("laminlabs/lamin-site-assets").get(1)
    ln.Artifact.connect("laminlabs/lamin-site-assets").get(uid="MqEaGU7fXvxNy61R0000")
    # cross-database query
    hemangioblast = bt.CellType.from_source(name="hemangioblast").save()
    artifact = (
        ln.Artifact.connect("laminlabs/lamin-dev")
        .filter(cell_types=hemangioblast)
        .first()
    )
    assert artifact is not None
    hemangioblast_dev = artifact.cell_types.get(name="hemangioblast")
    assert hemangioblast_dev.uid == hemangioblast.uid
    assert hemangioblast_dev.id != hemangioblast.id
    # query via list
    artifact_ref = (
        ln.Artifact.connect("laminlabs/lamin-dev")
        .filter(cell_types__in=[hemangioblast])
        .first()
    )
    assert artifact == artifact_ref
    # check that .using provided with the current intance does nothing
    assert ln.User.connect("lamindb-unit-tests-core").first()._state.db == "default"
    user = ln.setup.settings.user.handle
    assert (
        ln.User.connect(f"{user}/lamindb-unit-tests-core").first()._state.db
        == "default"
    )


def test_get_record_kwargs():
    assert _get_record_kwargs(ln.Feature) == [
        ("name", "str"),
        ("dtype", "DtypeStr | ULabel | Record | Registry | list[Registry] | FieldAttr"),
        ("type", "Feature | None"),
        ("is_type", "bool"),
        ("unit", "str | None"),
        ("description", "str | None"),
        ("synonyms", "str | None"),
        ("nullable", "bool | None"),
        (
            "default_value",
            "Any | None",
        ),
        ("coerce", "bool | None"),
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


def test_soft_delete_error():
    with pytest.raises(ValueError):
        ln.Storage.filter().first().delete(permanent=False)

    with pytest.raises(ValueError):
        ln.Branch.filter().first().delete(permanent=False)


def test_delete_return_value_permanent():
    """Test that permanent delete returns Django's natural return value."""
    # Test with ULabel (simple SQLRecord)
    ulabel = ln.ULabel(name="test-delete-return").save()
    result = ulabel.delete(permanent=True)
    assert isinstance(result, tuple)
    assert len(result) == 2
    deleted_count, deleted_dict = result
    assert deleted_count == 1
    assert isinstance(deleted_dict, dict)
    assert len(deleted_dict) > 0
    # Check that the registry name is in the dict
    # Django returns app_label.ClassName format
    registry_name = f"{ulabel._meta.app_label}.{ulabel.__class__.__name__}"
    assert registry_name in deleted_dict
    assert deleted_dict[registry_name] == 1


def test_unsaved_relationship_modification_attempts():
    af = ln.Artifact.from_dataframe(
        pd.DataFrame({"col1": [1, 2, 3], "col2": [4, 5, 6]}), description="testme"
    )

    new_label = ln.Record(name="testlabel").save()
    with pytest.raises(ValueError) as excinfo:
        af.records.add(new_label)

    assert (
        str(excinfo.value)
        == "You are trying to access the many-to-many relationships of an unsaved Artifact object. Please save it first using '.save()'."
    )

    new_label.delete(permanent=True)
    af.delete(permanent=True)


def test_failed_connect():
    with pytest.raises(ln.setup.errors.InstanceNotFoundError) as error:
        ln.Artifact.connect("laminlabs/lamindata-not-existing")
    assert error.exconly().startswith(
        "lamindb_setup.errors.InstanceNotFoundError: 'laminlabs/lamindata-not-existing' not found: 'instance-not-found'"
    )


def test_unsaved_model_different_instance():
    af = ln.Artifact.connect("laminlabs/lamindata").get(
        key="scrna/micro-macfarland2020.h5ad"
    )

    new_label = ln.Record(name="testlabel").save()
    with pytest.raises(ValueError) as excinfo:
        af.records.add(new_label)

    assert (
        str(excinfo.value)
        == "Cannot label a record from instance 'laminlabs/lamindata'. "
        "Please save the record first to your instance using '.save()'."
    )

    new_label.delete(permanent=True)
