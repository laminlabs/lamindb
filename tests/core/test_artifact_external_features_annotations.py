# ruff: noqa: F811

from datetime import date, datetime

import bionty as bt
import lamindb as ln
import pytest
from lamindb.errors import DoesNotExist, ValidationError
from lamindb.examples.datasets import mini_immuno
from lamindb.models._feature_manager import describe_features


# see test_record_basics.py for similar test for records
def test_artifact_features_add_remove_values():
    record_type1 = ln.Record(name="RecordType1", is_type=True).save()
    record_entity1 = ln.Record(name="entity1", type=record_type1).save()
    record_entity2 = ln.Record(name="entity2", type=record_type1).save()
    ulabel = ln.ULabel(name="test-ulabel").save()
    artifact = ln.Artifact(".gitignore", key="test-artifact").save()
    transform = ln.Transform(key="test-transform").save()
    run = ln.Run(transform, name="test-run").save()

    feature_str = ln.Feature(name="feature_str", dtype=str).save()
    feature_int = ln.Feature(name="feature_int", dtype=int).save()
    feature_datetime = ln.Feature(name="feature_datetime", dtype=datetime).save()
    feature_date = ln.Feature(name="feature_date", dtype=datetime.date).save()
    feature_dict = ln.Feature(name="feature_dict", dtype=dict).save()
    feature_type1 = ln.Feature(name="feature_type1", dtype=record_type1).save()
    feature_type1s = ln.Feature(name="feature_type1s", dtype=list[record_type1]).save()
    feature_ulabel = ln.Feature(name="feature_ulabel", dtype=ln.ULabel).save()
    feature_user = ln.Feature(name="feature_user", dtype=ln.User).save()
    feature_project = ln.Feature(name="feature_project", dtype=ln.Project).save()
    # feature_artifact = ln.Feature(name="feature_artifact", dtype=ln.Artifact).save()
    # feature_run = ln.Feature(name="feature_run", dtype=ln.Run.uid).save()
    feature_cell_line = ln.Feature(name="feature_cell_line", dtype=bt.CellLine).save()
    feature_cell_lines = ln.Feature(
        name="feature_cell_lines", dtype=list[bt.CellLine]
    ).save()
    feature_cl_ontology_id = ln.Feature(
        name="feature_cl_ontology_id", dtype=bt.CellLine.ontology_id
    ).save()

    test_artifact = ln.Artifact(".gitignore", key="test_artifact").save()
    test_project = ln.Project(name="test_project").save()
    hek293 = bt.CellLine.from_source(name="HEK293").save()
    a549 = bt.CellLine.from_source(name="A549 cell").save()

    # no schema validation

    test_values = {
        "feature_str": "a string value",
        "feature_int": 42,
        "feature_datetime": datetime(2024, 1, 1, 12, 0, 0),
        "feature_date": date(2024, 1, 1),
        "feature_dict": {"key": "value", "number": 123, "list": [1, 2, 3]},
        "feature_type1": "entity1",
        "feature_type1s": ["entity1", "entity2"],
        "feature_ulabel": "test-ulabel",
        "feature_user": ln.setup.settings.user.handle,
        "feature_project": "test_project",
        "feature_cell_line": "HEK293",
        "feature_cell_lines": ["HEK293", "A549 cell"],
        "feature_cl_ontology_id": "CLO:0001230",
        # "feature_artifact": "test-artifact",
        # "feature_run": run.uid,
    }

    test_artifact.features.add_values(test_values)
    assert test_artifact.features.get_values() == test_values

    # remove values

    test_artifact.features.remove_values("feature_int")
    test_values.pop("feature_int")
    assert test_artifact.features.get_values() == test_values

    test_artifact.features.remove_values("feature_date")
    test_values.pop("feature_date")
    assert test_artifact.features.get_values() == test_values

    test_artifact.features.remove_values("feature_type1")
    test_values.pop("feature_type1")
    assert test_artifact.features.get_values() == test_values

    test_artifact.features.remove_values("feature_type1s")
    test_values.pop("feature_type1s")
    assert test_artifact.features.get_values() == test_values

    test_artifact.features.remove_values("feature_ulabel")
    test_values.pop("feature_ulabel")
    assert test_artifact.features.get_values() == test_values

    # test passing a list

    test_artifact.features.remove_values(["feature_cell_line", "feature_user"])
    test_values.pop("feature_cell_line")
    test_values.pop("feature_user")
    assert test_artifact.features.get_values() == test_values

    # test_artifact.features.remove_values("feature_artifact")
    # test_values.pop("feature_artifact")
    # assert test_artifact.features.get_values() == test_values

    # test_artifact.features.remove_values("feature_run")
    # test_values.pop("feature_run")
    # assert test_artifact.features.get_values() == test_values

    # test passing None has no effect, does not lead to annotation

    test_artifact.features.add_values({"feature_int": None, "feature_type1": None})
    assert test_artifact.features.get_values() == test_values

    # test bulk removal

    assert list(test_values.keys()) == [
        "feature_str",
        "feature_datetime",
        "feature_dict",
        "feature_project",
        "feature_cell_lines",
        "feature_cl_ontology_id",
    ]
    test_artifact.features.remove_values()
    test_values = {}
    assert test_artifact.features.get_values() == test_values

    # test passing ISO-format date string for date

    test_artifact.features.add_values({"feature_date": "2024-01-01"})
    test_values["feature_date"] = date(2024, 1, 1)
    assert test_artifact.features.get_values() == test_values

    # test add_values() when there is already something there

    test_artifact.features.add_values({"feature_date": "2024-02-01"})
    test_values["feature_date"] = {date(2024, 1, 1), date(2024, 2, 1)}
    test_artifact.features.add_values({"feature_str": "a string value"})
    test_values["feature_str"] = "a string value"
    assert test_artifact.features.get_values() == test_values

    # test set_values()

    test_values = {}
    test_values["feature_date"] = date(2024, 3, 1)
    test_artifact.features.set_values({"feature_date": "2024-03-01"})
    assert test_artifact.features.get_values() == test_values

    # schema validation

    feature_str = ln.Feature.get(name="feature_str")
    feature_int = ln.Feature.get(name="feature_int")
    schema = ln.Schema([feature_str, feature_int], name="test_schema").save()
    with pytest.raises(ln.errors.ValidationError) as error:
        test_artifact.features.add_values({"feature_type1": "entity1"}, schema=schema)
    assert "COLUMN_NOT_IN_DATAFRAME" in error.exconly()
    schema.delete(permanent=True)

    # test with list of strings

    schema = ln.Schema([feature_cell_lines], name="test_schema2").save()
    test_artifact.features.add_values(
        {"feature_cell_lines": ["HEK293", "A549 cell"]}, schema=schema
    )
    schema.delete(permanent=True)

    # test with list of records (rather than passing strings)

    schema = ln.Schema([feature_cell_lines], name="test_schema2").save()
    test_artifact.features.add_values(
        {"feature_cell_lines": [a549, hek293]}, schema=schema
    )
    schema.delete(permanent=True)

    # clean up rest

    test_artifact.delete(permanent=True)
    feature_str.delete(permanent=True)
    feature_int.delete(permanent=True)
    feature_datetime.delete(permanent=True)
    feature_date.delete(permanent=True)
    feature_type1.delete(permanent=True)
    feature_type1s.delete(permanent=True)
    feature_user.delete(permanent=True)
    feature_project.delete(permanent=True)
    feature_dict.delete(permanent=True)
    # feature_artifact.delete(permanent=True)
    # feature_run.delete(permanent=True)
    feature_ulabel.delete(permanent=True)
    feature_cell_lines.delete(permanent=True)
    record_entity1.delete(permanent=True)
    record_entity2.delete(permanent=True)
    record_type1.delete(permanent=True)
    test_project.delete(permanent=True)
    feature_cell_line.delete(permanent=True)
    feature_cl_ontology_id.delete(permanent=True)
    hek293.delete(permanent=True)
    a549.delete(permanent=True)
    ulabel.delete(permanent=True)
    artifact.delete(permanent=True)
    run.delete(permanent=True)
    transform.delete(permanent=True)


def test_features_add_with_schema():
    df = mini_immuno.get_dataset1(otype="DataFrame")
    artifact = ln.Artifact.from_dataframe(df, description="test dataset").save()

    species = ln.Feature(name="species", dtype="str").save()
    split = ln.Feature(name="split", dtype="str").save()
    schema = ln.Schema([species, split]).save()

    with pytest.raises(ValidationError) as e:
        artifact.features.add_values({"doesnot": "exist"}, schema=schema)
    assert "column 'split' not in dataframe" in str(e.value)

    artifact.features.add_values({"species": "bird", "split": "train"}, schema=schema)
    artifact.save()

    assert artifact.features.get_values() == {"species": "bird", "split": "train"}

    artifact.delete(permanent=True)
    schema.delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


def test_features_add_remove_error_behavior():
    adata = ln.examples.datasets.anndata_with_obs()
    artifact = ln.Artifact.from_anndata(adata, description="test").save()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"experiment": "Experiment 1"})
    assert (
        error.exconly()
        == """lamindb.errors.ValidationError: These keys could not be validated: ['experiment']
Here is how to create a feature:

  ln.Feature(name='experiment', dtype='cat ? str').save()"""
    )
    ln.Feature(name="experiment", dtype=ln.Record).save()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"experiment": "Experiment 1"})
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: These values could not be validated:"
    )
    ln.Record(name="Experiment 1").save()
    # now add the label with the feature and make sure that it has the feature annotation
    artifact.features.add_values({"experiment": "Experiment 1"})
    assert artifact.links_record.get().record.name == "Experiment 1"
    assert artifact.links_record.get().feature.name == "experiment"
    # repeat
    artifact.features.add_values({"experiment": "Experiment 1"})
    assert artifact.links_record.get().record.name == "Experiment 1"

    # numerical feature
    temperature = ln.Feature(name="temperature", dtype="cat").save()
    with pytest.raises(TypeError) as error:
        artifact.features.add_values({"temperature": 27.2})
    assert error.exconly().startswith(
        "TypeError: Value for feature 'temperature' with dtype 'cat' must be a string or record"
    )
    temperature.dtype = "num"
    temperature.save()
    artifact.features.add_values({"temperature": 27.2})
    assert artifact._feature_values.first().value == 27.2

    # datetime feature
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"date_of_experiment": "2024-12-01"})
    assert (
        error.exconly()
        == """lamindb.errors.ValidationError: These keys could not be validated: ['date_of_experiment']
Here is how to create a feature:

  ln.Feature(name='date_of_experiment', dtype='date').save()"""
    )

    ln.Feature(name="date_of_experiment", dtype="date").save()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"date_of_experiment": "Typo2024-12-01"})
    assert (
        error.exconly()
        == """lamindb.errors.ValidationError: Expected dtype for 'date_of_experiment' is 'date', got 'cat ? str'"""
    )
    artifact.features.add_values({"date_of_experiment": "2024-12-01"})

    ln.Feature(name="datetime_of_experiment", dtype="datetime").save()
    artifact.features.add_values({"datetime_of_experiment": "2024-12-01 00:00:00"})

    # bionty feature
    mouse = bt.Organism.from_source(name="mouse")
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"organism": mouse})
    assert (
        error.exconly()
        == """lamindb.errors.ValidationError: These keys could not be validated: ['organism']
Here is how to create a feature:

  ln.Feature(name='organism', dtype='cat[bionty.Organism]').save()"""
    )
    ln.Feature(name="organism", dtype=bt.Organism).save()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"organism": mouse})
    assert (
        # ensure the label is saved
        error.exconly().startswith("lamindb.errors.ValidationError: Please save")
    )
    mouse.save()
    artifact.features.add_values({"organism": mouse})
    assert artifact.organisms.get().name == "mouse"

    # lists of records
    diseases = bt.Disease.from_values(
        ["MONDO:0004975", "MONDO:0004980"], field=bt.Disease.ontology_id
    )
    ln.save(diseases)
    ln.Feature(name="disease", dtype="cat[bionty.Disease.ontology_id]").save()
    artifact.features.add_values({"disease": diseases})
    assert len(artifact.diseases.filter()) == 2
    # check get_values returns ontology_ids as specified in the feature dtype
    assert artifact.features.get_values()["disease"] == {
        "MONDO:0004975",
        "MONDO:0004980",
    }

    # big dictionary of everything
    features = {
        "experiment": [  # we're testing iterable annotation here
            "Experiment 2",
            "Experiment 1",
        ],
        "project": "project_1",
        "is_validated": True,
        "cell_type_by_expert": "T cell",
        "temperature": 100.0,
        "donor": "U0123",
    }
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values(features)
    assert (
        error.exconly()
        == """\
lamindb.errors.ValidationError: These keys could not be validated: ['project', 'is_validated', 'cell_type_by_expert', 'donor']
Here is how to create a feature:

  ln.Feature(name='project', dtype='cat ? str').save()
  ln.Feature(name='is_validated', dtype='bool').save()
  ln.Feature(name='cell_type_by_expert', dtype='cat ? str').save()
  ln.Feature(name='donor', dtype='cat ? str').save()"""
    )

    ln.Feature(name="project", dtype=ln.Record).save()
    ln.Feature(name="is_validated", dtype=bool).save()
    ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save()
    ln.Feature(name="donor", dtype=ln.Record).save()

    with pytest.raises(ValidationError) as error:
        artifact.features.add_values(features)
        error_msg = error.exconly()

        assert (
            "lamindb.errors.ValidationError: These values could not be validated:"
            in error_msg
        )
        assert "Here is how to create records for them:" in error_msg

        expected_values = {
            "Record": ["project_1", "U0123", "Experiment 2"],
            "bionty.CellType": ["T cell"],
        }

        for key, values in expected_values.items():
            assert f"'{key}':" in error_msg
            for value in values:
                assert value in error_msg
            assert f"{key.split('.')[-1]}.from_values(" in error_msg

        assert "create=True).save()" in error_msg

    ln.Record.from_values(["Experiment 2", "project_1", "U0123"], create=True).save()
    bt.CellType.from_source(name="T cell").save()
    print("validate", bt.CellType.validate(["T cell"]))

    artifact.features.add_values(features)
    assert set(artifact._feature_values.all().values_list("value", flat=True)) == {
        27.2,
        True,
        100.0,
        "2024-12-01",
        "2024-12-01T00:00:00",
    }

    assert ln.Artifact.get(_feature_values__value=27.2)

    assert artifact.features.get_values() == {
        "disease": {"MONDO:0004975", "MONDO:0004980"},
        "experiment": {"Experiment 1", "Experiment 2"},
        "project": "project_1",
        "cell_type_by_expert": "T cell",
        "donor": "U0123",
        "organism": "mouse",
        "is_validated": True,
        "temperature": {27.2, 100.0},
        "date_of_experiment": date(2024, 12, 1),
        "datetime_of_experiment": datetime(2024, 12, 1, 0, 0, 0),
    }
    # hard to test because of italic formatting
    _, external_features_tree = describe_features(artifact)
    assert external_features_tree.label.plain == "Features"
    assert len(external_features_tree.children[0].label.columns) == 3
    assert len(external_features_tree.children[0].label.rows) == 10
    assert external_features_tree.children[0].label.columns[0]._cells == [
        "cell_type_by_expert",
        "disease",
        "donor",
        "experiment",
        "organism",
        "project",
        "date_of_experiment",
        "datetime_of_experiment",
        "is_validated",
        "temperature",
    ]
    dtypes_display = [
        i.plain for i in external_features_tree.children[0].label.columns[1]._cells
    ]
    assert dtypes_display == [
        "bionty.CellType",
        "bionty.Disease.ontology_id",
        "Record",
        "Record",
        "bionty.Organism",
        "Record",
        "date",
        "datetime",
        "bool",
        "num",
    ]
    assert external_features_tree.children[0].label.columns[2]._cells == [
        "T cell",
        "MONDO:0004975, MONDO:0004980",
        "U0123",
        "Experiment 1, Experiment 2",
        "mouse",
        "project_1",
        "2024-12-01",
        "2024-12-01 00:00:00",
        "True",
        "27.2, 100.0",
    ]
    # repeat
    artifact.features.add_values(features)
    assert set(artifact._feature_values.all().values_list("value", flat=True)) == {
        27.2,
        True,
        100.0,
        "2024-12-01",
        "2024-12-01T00:00:00",
    }

    with pytest.raises(ln.errors.InvalidArgument) as error:
        ln.Artifact.filter(temperature_with_typo=100.0, project="project_1").one()
    assert error.exconly().startswith(
        "lamindb.errors.InvalidArgument: You can query either by available fields:"
    )

    ln.Artifact.filter(temperature=100.0)
    ln.Artifact.filter(project="project_1")
    ln.Artifact.filter(is_validated=True)
    ln.Artifact.filter(temperature=100.0, project="project_1", donor="U0123").one()
    # for bionty
    assert (
        artifact == ln.Artifact.filter(disease=diseases[0]).one()
    )  # value is a record
    assert (
        artifact == ln.Artifact.filter(disease="MONDO:0004975").one()
    )  # value is a string
    assert artifact == ln.Artifact.filter(disease__contains="0004975").one()

    # test not finding the Record
    with pytest.raises(DoesNotExist) as error:
        ln.Artifact.filter(project="project__1")
    assert error.exconly().startswith(
        "lamindb.errors.DoesNotExist: Did not find a Record matching"
    )

    # test comparator
    assert artifact == ln.Artifact.filter(experiment__contains="ment 1").one()
    # due to the __in comparator, we get the same artifact twice below
    # print(ln.Artifact.to_dataframe(features=["experiment"]))
    # print(ln.Artifact.filter(experiment__contains="Experi").to_dataframe(features=["experiment"]))
    assert len(ln.Artifact.filter(experiment__contains="Experi")) == 2
    assert ln.Artifact.filter(temperature__lt=21).one_or_none() is None
    assert len(ln.Artifact.filter(temperature__gt=21)) >= 1

    # test remove_values
    artifact.features.remove_values("date_of_experiment")
    alzheimer = bt.Disease.get(name="Alzheimer disease")
    artifact.features.remove_values("disease", value=alzheimer)
    values = artifact.features.get_values()
    assert "date_of_experiment" not in values
    assert "MONDO:0004975" not in values["disease"]

    # test annotate with dictionaries multiple times
    ln.Feature(name="study_metadata", dtype=dict).save()
    artifact.features.add_values({"study_metadata": {"detail1": "123", "detail2": 1}})

    # delete everything we created
    artifact.delete(permanent=True)
    ln.Record.filter().delete(permanent=True)
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    bt.Gene.filter().delete(permanent=True)
    bt.Organism.filter().delete(permanent=True)
    bt.Disease.filter().delete(permanent=True)


def test_add_remove_list_features(ccaplog):
    feature = ln.Feature(name="list_of_str", dtype=list[str]).save()
    artifact = ln.Artifact(".gitignore", key=".gitignore").save()
    artifact.features.add_values({"list_of_str": ["1", "2", "3"]})
    assert artifact.features.get_values() == {"list_of_str": ["1", "2", "3"]}
    # remove a non-linked value, this should do nothing but print a warning
    artifact.features.remove_values("list_of_str", value="4")
    assert "no feature 'list_of_str' with value '4' found" in ccaplog.text
    # list of categories feature
    cell_types_feature = ln.Feature(
        name="cell_types", dtype="list[cat[bionty.CellType]]"
    ).save()
    bt.CellType.from_values(["T cell", "B cell"]).save()
    artifact.features.add_values({"cell_types": ["T cell", "B cell"]})
    assert set(artifact.features.get_values()["cell_types"]) == {"B cell", "T cell"}
    # passing value works here because we are linking each of the cell types in the list individually
    # in comparison to passing a list of numbers above
    t_cell = bt.CellType.get(name="T cell")
    artifact.features.remove_values("cell_types", value=t_cell)
    assert artifact.features.get_values()["cell_types"] == ["B cell"]
    # remove a non-linked value, this should print a warning but do nothing
    artifact.features.remove_values("cell_types", value=t_cell.parents.first())
    assert "no feature 'cell_types' with value CellType(" in ccaplog.text
    # remove the entire linked feature
    artifact.features.remove_values("cell_types")
    assert "cell_types" not in artifact.features.get_values()

    # clean up
    artifact.delete(permanent=True)
    assert ln.models.FeatureValue.filter(feature__name="list_of_str").count() == 1
    feature.delete(permanent=True)
    assert ln.models.FeatureValue.filter(feature__name="list_of_str").count() == 0
    cell_types_feature.delete(permanent=True)
    bt.CellType.filter().delete(permanent=True)


def test_add_list_of_cat_features():
    type_1 = ln.Record(name="Type 1", is_type=True).save()
    for label in ["label 1", "label 2", "label 3"]:
        ln.Record(name=label, type=type_1).save()
    feat1 = ln.Feature(
        name="single_label_of_type1", dtype=type_1, nullable=False
    ).save()
    feat2 = ln.Feature(
        name="list_of_labels_of_type1", dtype=list[type_1], nullable=False
    ).save()
    schema = ln.Schema(name="Test schema", features=[feat1, feat2]).save()
    artifact = ln.Artifact(
        ".gitignore",
        key=".gitignore",
    ).save()
    # now just use add_values()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values(
            {
                "single_label_of_type1": "invalid",
            }
        )
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: These values could not be validated: {'Record': ('name', ['invalid'])}"
    )
    # now for list of labels
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values(
            {
                "list_of_labels_of_type1": ["invalid", "invalid2"],
            }
        )
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: These values could not be validated: {'Record': ('name', ['invalid', 'invalid2'])}"
    )
    artifact.delete(permanent=True)
    # now with schema
    artifact = ln.Artifact(
        ".gitignore",
        key=".gitignore",
        schema=schema,
        features={
            "single_label_of_type1": "label 1",
            "list_of_labels_of_type1": ["label 1", "label 2"],
        },
    ).save()
    with pytest.raises(ValueError) as error:
        artifact.features.add_values(
            {
                "single_label_of_type1": "invalid",
            }
        )
    assert "Cannot add values if artifact has external schema." in error.exconly()

    artifact.delete(permanent=True)
    schema.delete(permanent=True)
    feat1.delete(permanent=True)
    feat2.delete(permanent=True)
    type_1.records.all().delete(permanent=True)
    type_1.delete(permanent=True)
