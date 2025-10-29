# ruff: noqa: F811

import datetime

import bionty as bt
import lamindb as ln
import pytest
from lamindb.errors import DoesNotExist, ValidationError
from lamindb.examples.datasets import mini_immuno
from lamindb.models._feature_manager import describe_features


@pytest.fixture(scope="module")
def adata():
    adata = ln.examples.datasets.anndata_with_obs()
    # add another column
    adata.obs["cell_type_by_expert"] = adata.obs["cell_type"]
    adata.obs.loc["obs0", "cell_type_by_expert"] = "B cell"
    return adata


def test_features_add():
    df = mini_immuno.get_dataset1(otype="DataFrame")
    artifact = ln.Artifact.from_dataframe(df, description="test dataset").save()
    with pytest.raises(ValidationError) as err:
        artifact.features.add_values({"perturbation": df.perturbation.unique()})
    assert (
        err.exconly()
        == """lamindb.errors.ValidationError: These keys could not be validated: ['perturbation']
Here is how to create a feature:

  ln.Feature(name='perturbation', dtype='cat').save()"""
    )

    perturbation_feature = ln.Feature(name="perturbation", dtype=ln.Record).save()
    records = ln.Record.from_values(["DMSO", "IFNG"], create=True).save()
    artifact.features.add_values({"perturbation": df.perturbation.unique()})
    assert artifact in ln.Artifact.filter(perturbation__isnull=False)
    assert artifact not in ln.Artifact.filter(perturbation__isnull=True)

    # list of bionty features
    organisms_feature = ln.Feature(name="organisms", dtype=list[bt.Organism]).save()
    mouse = bt.Organism.from_source(name="mouse").save()
    artifact.features.add_values({"organisms": [mouse]})
    assert artifact.features.get_values()["organisms"] == ["mouse"]

    artifact.delete(permanent=True)
    organisms_feature.delete(permanent=True)
    mouse.delete(permanent=True)
    perturbation_feature.delete(permanent=True)
    for record in records:
        record.delete(permanent=True)


def test_features_add_external():
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


# below the test for annotating with feature values
def test_features_add_remove(adata, ccaplog):
    artifact = ln.Artifact.from_anndata(adata, description="test").save()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"experiment": "Experiment 1"})
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: These keys could not be validated:"
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
        "date_of_experiment": datetime.date(2024, 12, 1),
        "datetime_of_experiment": datetime.datetime(2024, 12, 1, 0, 0, 0),
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
    # first end-to-end test of adding labels and passing a schema
    # this calls features.add_values() under-the-hood
    artifact = ln.Artifact(
        ".gitignore",
        key=".gitignore",
        schema=schema,
        features={
            "single_label_of_type1": "label 1",
            "list_of_labels_of_type1": ["label 1", "label 2"],
        },
    ).save()
    # now just use add_values()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values(
            {
                "single_label_of_type1": "invalid",
            }
        )
    print(error.exconly())
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
    print(error.exconly())
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: These values could not be validated: {'Record': ('name', ['invalid', 'invalid2'])}"
    )

    artifact.delete(permanent=True)
    schema.delete(permanent=True)
    feat1.delete(permanent=True)
    feat2.delete(permanent=True)
    type_1.records.all().delete(permanent=True)
    type_1.delete(permanent=True)
