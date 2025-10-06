# ruff: noqa: F811

import datetime
from pathlib import Path

import bionty as bt
import lamindb as ln
import pytest
from _dataset_fixtures import (  # noqa
    get_mini_csv,
)
from lamindb.errors import DoesNotExist, ValidationError
from lamindb.examples.datasets import small_dataset1
from lamindb.models._feature_manager import describe_features
from lamindb.models._label_manager import format_rich_tree
from lamindb.models.artifact import add_labels


@pytest.fixture(scope="module")
def adata():
    adata = ln.core.datasets.anndata_with_obs()
    # add another column
    adata.obs["cell_type_by_expert"] = adata.obs["cell_type"]
    adata.obs.loc["obs0", "cell_type_by_expert"] = "B cell"
    return adata


def test_features_add():
    df = small_dataset1(otype="DataFrame")
    artifact = ln.Artifact.from_dataframe(df, description="test dataset").save()
    with pytest.raises(ValidationError) as err:
        artifact.features.add_values({"perturbation": df.perturbation.unique()})
    assert (
        err.exconly()
        == """lamindb.errors.ValidationError: These keys could not be validated: ['perturbation']
Here is how to create a feature:

  ln.Feature(name='perturbation', dtype='cat').save()"""
    )

    ln.Feature(name="perturbation", dtype=ln.Record).save()
    ln.Record.from_values(["DMSO", "IFNG"], create=True).save()
    artifact.features.add_values({"perturbation": df.perturbation.unique()})
    assert artifact in ln.Artifact.filter(perturbation__isnull=False)
    assert artifact not in ln.Artifact.filter(perturbation__isnull=True)

    artifact.delete(permanent=True)
    ln.Record.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


def test_features_add_external():
    df = small_dataset1(otype="DataFrame")
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
def test_features_add_remove(adata):
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
    experiment_label = ln.Record(name="Experiment 1").save()
    # add the label without the feature first
    artifact.records.add(experiment_label)
    assert artifact.links_record.get().record.name == "Experiment 1"
    assert artifact.links_record.get().feature is None

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
    ln.Feature(name="organism", dtype="cat[bionty.Organism]").save()
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
    ln.Feature(name="disease", dtype="cat[bionty.Disease]").save()
    artifact.features.add_values({"disease": diseases})
    assert len(artifact.diseases.filter().all()) == 2

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
        "disease": {"Alzheimer disease", "atopic eczema"},
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
    tree = describe_features(artifact)
    format_rich_tree(tree)
    assert tree.children[0].label.plain == "Linked features"
    assert len(tree.children[0].children[0].label.columns) == 3
    assert len(tree.children[0].children[0].label.rows) == 10
    assert tree.children[0].children[0].label.columns[0]._cells == [
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
    types = [i.plain for i in tree.children[0].children[0].label.columns[1]._cells]
    assert types == [
        "cat[bionty.CellType]",
        "cat[bionty.Disease]",
        "cat[Record]",
        "cat[Record]",
        "cat[bionty.Organism]",
        "cat[Record]",
        "date",
        "datetime",
        "bool",
        "num",
    ]
    assert tree.children[0].children[0].label.columns[2]._cells == [
        "T cell",
        "Alzheimer disease, atopic eczema",
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
        artifact == ln.Artifact.filter(disease="Alzheimer disease").one()
    )  # value is a string
    assert artifact == ln.Artifact.filter(disease__contains="Alzheimer").one()

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
    assert len(ln.Artifact.filter(experiment__contains="Experi").all()) == 2
    assert ln.Artifact.filter(temperature__lt=21).one_or_none() is None
    assert len(ln.Artifact.filter(temperature__gt=21).all()) >= 1

    # test remove_values
    artifact.features.remove_values("date_of_experiment")
    alzheimer = bt.Disease.get(name="Alzheimer disease")
    artifact.features.remove_values("disease", value=alzheimer)
    values = artifact.features.get_values()
    assert "date_of_experiment" not in values
    assert "Alzheimer disease" not in values["disease"]

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


# most underlying logic here is comprehensively tested in test_context
def test_params_add():
    path = Path("mymodel.pt")
    path.touch()
    artifact = ln.Artifact("mymodel.pt", kind="model", description="hello").save()
    ln.Feature(name="learning_rate", dtype="float").save()
    ln.Feature(name="quantification", dtype="dict").save()
    artifact.features.add_values({"learning_rate": 0.01})
    artifact.features.add_values(
        {
            "quantification": {
                "name": "mcquant",
                "container": "labsyspharm/quantification",
            }
        }
    )
    assert artifact.features.get_values() == {
        "learning_rate": 0.01,
        "quantification": {
            "name": "mcquant",
            "container": "labsyspharm/quantification",
        },
    }
    # test describe params
    tree = describe_features(artifact)
    assert tree.label.plain == "Artifact .pt"
    assert tree.children[0].label.plain == "Linked features"
    assert len(tree.children[0].children[0].label.columns) == 3
    assert tree.children[0].children[0].label.columns[0]._cells == [
        "learning_rate",
        "quantification",
    ]
    assert tree.children[0].children[0].label.columns[1]._cells[0].plain == "float"
    assert tree.children[0].children[0].label.columns[1]._cells[1].plain == "dict"
    assert tree.children[0].children[0].label.columns[2]._cells == [
        "0.01",
        "{'name': 'mcquant', 'container': 'labsyspharm/quantification'}",
    ]
    artifact.describe()
    artifact.delete(permanent=True)
    path.unlink()


def test_labels_add(adata):
    label = ln.Record(name="Experiment 1")
    artifact = ln.Artifact.from_anndata(adata, description="test").save()
    experiment = ln.Feature(name="experiment", dtype=ln.Record)
    with pytest.raises(ValueError) as error:
        artifact.labels.add("experiment_1", experiment)
    assert (
        error.exconly()
        == "ValueError: Please pass a record (a `SQLRecord` object), not a string, e.g.,"
        " via: label = ln.Record(name='experiment_1')"
    )
    with pytest.raises(ValidationError) as error:
        artifact.labels.add(label, experiment)
    assert "not validated. If it looks correct: record.save()" in error.exconly()
    label.save()
    with pytest.raises(TypeError) as error:
        artifact.labels.add(label, "experiment 1")
    with pytest.raises(ValidationError) as error:
        artifact.labels.add(label, feature=experiment)
    assert (
        error.exconly()
        == "lamindb.errors.ValidationError: Feature not validated. If it looks"
        " correct: ln.Feature(name='experiment', type='cat[Record]').save()"
    )
    experiment.save()

    # try to pass list of length zero
    artifact.labels.add([], feature=experiment)
    # now pass a single label
    artifact.labels.add(label, feature=experiment)
    # check that the feature was updated with type = "Record"
    feature = ln.Feature.get(name="experiment")
    assert feature.dtype == "cat[Record]"
    with pytest.raises(TypeError):
        experiments = artifact.labels.get("experiment")
    # check that the label is there, it's exactly one label with name "Experiment 1"
    experiments = artifact.labels.get(experiment)
    assert experiments.one().name == "Experiment 1"

    # try adding the same label again, nothing should happen
    artifact.labels.add(label, feature=experiment)
    # check that the label is there, it's exactly one label with name "Experiment 1"
    experiments = artifact.labels.get(experiment)
    assert experiments.get().name == "Experiment 1"

    # running from_values to load validated label records under the hood
    experiment = ln.Feature(name="experiment_with_reg", dtype="cat[Record]").save()
    ln.Record(name="Experiment 2").save()
    artifact.labels.add("Experiment 2", experiment)
    experiments = artifact.labels.get(experiment)
    assert experiments.get().name == "Experiment 2"

    # now, try adding a new label
    project = ln.Record(name="project 1").save()
    ln.Feature(name="project", dtype=ln.Record).save()
    features = ln.Feature.lookup()
    artifact.labels.add(project, feature=features.project)
    # check that the label is there, it's exactly one label with name "Experiment 1"
    projects = artifact.labels.get(features.project)
    assert projects.get().name == "project 1"

    # test add_from
    adata2 = adata.copy()
    adata2.uns["mutated"] = True
    artifact2 = ln.Artifact(adata2, description="My new artifact").save()

    artifact2.labels.add_from(artifact)
    experiments = artifact2.labels.get(experiment)
    assert experiments.get().name == "Experiment 2"

    artifact2.delete(permanent=True)
    artifact.delete(permanent=True)
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    ln.Record.filter().delete(permanent=True)


def test_add_labels_using_anndata(adata):
    organism = bt.Organism.from_source(name="mouse")
    cell_types = [bt.CellType(name=name) for name in adata.obs["cell_type"].unique()]
    ln.save(cell_types)
    inspector = bt.CellType.inspect(adata.obs["cell_type_by_expert"].unique())
    ln.save([bt.CellType(name=name) for name in inspector.non_validated])
    cell_types_from_expert = bt.CellType.from_values(
        adata.obs["cell_type_by_expert"].unique()
    )
    actual_tissues = [bt.Tissue(name=name) for name in adata.obs["tissue"].unique()]
    organoid = ln.Record(name="organoid")
    tissues = actual_tissues + [organoid]
    ln.save(tissues)

    bt.settings.organism = "human"

    # clean up DB state
    organism_feature = ln.Feature.filter(name="organism").one_or_none()
    if organism_feature is not None:
        organism_feature.delete(permanent=True)
    artifact = ln.Artifact.filter(description="Mini adata").one_or_none()
    if artifact is not None:
        artifact.delete(permanent=True, storage=True)
    ln.Schema.filter().delete(permanent=True)

    # try to construct without registering metadata features
    artifact = ln.Artifact.from_anndata(adata, description="Mini adata")
    if not artifact._state.adding:
        artifact.delete(permanent=True)  # make sure we get a fresh one
        artifact = ln.Artifact.from_anndata(adata, description="Mini adata")
    # add feature set without saving file
    feature_name_feature = ln.Feature(name="feature name", dtype="cat[Record]").save()
    schema = ln.Schema(features=[feature_name_feature])
    with pytest.raises(ValueError) as error:
        artifact.features._add_schema(schema, slot="random")
    assert (
        error.exconly()
        == "ValueError: Please save the artifact or collection before adding a feature"
        " set!"
    )

    # now register features we want to validate
    # (we are not interested in cell_type_id, here)
    ln.save(
        ln.Feature.from_dataframe(
            adata.obs[["cell_type", "tissue", "cell_type_by_expert", "disease"]]
        )
    )
    artifact = ln.Artifact.from_anndata(adata, description="Mini adata")
    ln.Feature(name="organism", dtype="cat[bionty.Organism]").save()
    features = ln.Feature.lookup()
    with pytest.raises(ValueError) as error:
        artifact.labels.add(organism, feature=features.organism)
    assert (
        error.exconly()
        == "ValueError: Please save the artifact/collection before adding a label!"
    )
    artifact.save()

    # link features
    artifact.features._add_set_from_anndata(var_field=bt.Gene.ensembl_gene_id)

    # check the basic construction of the feature set based on obs
    schema_obs = artifact.feature_sets.filter(
        itype="Feature", _links_artifact__slot="obs"
    ).one()
    assert schema_obs.n == 4
    assert "organism" not in schema_obs.features.to_list("name")

    # now, we add organism and run checks
    features = ln.Feature.lookup()
    with pytest.raises(ln.errors.ValidationError):
        artifact.labels.add(organism, feature=features.organism)
    organism.save()
    artifact.labels.add(organism, feature=features.organism)
    organism_link = artifact.links_organism.first()
    assert organism_link.organism.name == "mouse"
    assert organism_link.feature.name == "organism"
    feature = ln.Feature.get(name="organism")
    assert feature.dtype == "cat[bionty.Organism]"
    schema_obs = artifact.feature_sets.filter(
        itype="Feature", _links_artifact__slot="obs"
    ).one()
    assert schema_obs.n == 4

    # now we add cell types & tissues and run checks
    ln.Feature(name="cell_type", dtype="cat").save()
    ln.Feature(name="cell_type_by_expert", dtype="cat").save()
    ln.Feature(name="tissue", dtype="cat").save()
    add_labels(artifact, cell_types, feature=features.cell_type, from_curator=True)
    add_labels(
        artifact,
        cell_types_from_expert,
        feature=features.cell_type_by_expert,
        from_curator=True,
    )
    with pytest.raises(ValidationError) as err:
        add_labels(artifact, tissues, feature=features.tissue, from_curator=True)
    assert (
        err.exconly()
        == "lamindb.errors.ValidationError: Label type Record is not valid for Feature(name='tissue', dtype='cat[bionty.Tissue]'), consider updating to dtype='cat[bionty.Tissue|Record]'"
    )
    tissue = ln.Feature.get(name="tissue")
    tissue.dtype = "cat[bionty.Tissue|Record]"
    tissue.save()
    add_labels(artifact, tissues, feature=tissue, from_curator=True)
    feature = ln.Feature.get(name="cell_type")
    assert feature.dtype == "cat[bionty.CellType]"
    feature = ln.Feature.get(name="cell_type_by_expert")
    assert feature.dtype == "cat[bionty.CellType]"
    feature = ln.Feature.get(name="tissue")
    assert feature.dtype == "cat[bionty.Tissue|Record]"
    diseases = [ln.Record(name=name) for name in adata.obs["disease"].unique()]
    ln.save(diseases)
    add_labels(artifact, diseases, feature=features.disease, from_curator=True)
    df = artifact.features.slots["obs"].features.to_dataframe()
    assert set(df["name"]) == {
        "cell_type",
        "disease",
        "tissue",
        "cell_type_by_expert",
    }
    assert set(df["dtype"]) == {
        "cat[bionty.CellType]",
        "cat[Record]",
        "cat[bionty.Tissue|Record]",
    }

    # now, let's add another feature to ext
    experiment_1 = ln.Record(name="experiment_1").save()
    ln.Feature(name="experiment", dtype="cat").save()
    features = ln.Feature.lookup()
    artifact.labels.add(experiment_1, feature=features.experiment)
    # TODO: replace the following with an updated test
    # df = artifact.features["external"].to_dataframe()
    # assert set(df["name"]) == {
    #     "organism",
    #     "experiment",
    # }
    # assert set(df["dtype"]) == {"cat[bionty.Organism]", "cat[Record]"}

    assert set(artifact.labels.get(features.experiment).to_list("name")) == {
        "experiment_1"
    }
    assert set(artifact.labels.get(features.disease).to_list("name")) == {
        "chronic kidney disease",
        "Alzheimer disease",
        "liver lymphoma",
        "cardiac ventricle disorder",
    }
    assert set(artifact.labels.get(features.organism).to_list("name")) == {"mouse"}
    assert set(
        artifact.labels.get(features.tissue)["bionty.Tissue"].to_list("name")
    ) == {
        "liver",
        "heart",
        "kidney",
        "brain",
    }
    assert set(artifact.labels.get(features.tissue)["Record"].to_list("name")) == {
        "organoid",
    }
    # currently, we can't stratify the two cases below
    assert set(artifact.labels.get(features.cell_type).to_list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(artifact.labels.get(features.cell_type, flat_names=True)) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(artifact.labels.get(features.cell_type_by_expert).to_list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert experiment_1 in artifact.records.all()

    # call describe
    artifact.describe()

    # clean up
    artifact.delete(permanent=True)
    bt.Gene.filter().delete(permanent=True)
    bt.Organism.filter().delete(permanent=True)
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    bt.CellType.filter().delete(permanent=True)
    bt.Tissue.filter().delete(permanent=True)
    bt.Disease.filter().delete(permanent=True)
    ln.Record.filter().delete(permanent=True)


def test_labels_get(get_mini_csv):
    artifact = ln.Artifact(get_mini_csv, description="test")
    # feature doesn't exist
    with pytest.raises(TypeError):
        artifact.labels.get("x")
    # no linked labels
    feature_name_feature = ln.Feature(name="feature name", dtype="cat")
    feature_name_feature.save()
    schema = ln.Schema(features=[feature_name_feature])
    schema.save()
    artifact.save()
    assert str(artifact.features) == "no linked features"
    # test for deprecated add_schema
    artifact.features._add_schema(schema, slot="random")
    assert artifact.feature_sets.first() == schema
    artifact.delete(permanent=True, storage=True)
    schema.delete(permanent=True)
    feature_name_feature.delete(permanent=True)


@pytest.fixture
def get_test_artifacts():
    with open("./default_storage_unit_core/test-inherit1", "w") as f:
        f.write("artifact1")
    with open("./default_storage_unit_core/test-inherit2", "w") as f:
        f.write("artifact2")
    artifact1 = ln.Artifact("./default_storage_unit_core/test-inherit1")
    artifact1.save()
    artifact2 = ln.Artifact("./default_storage_unit_core/test-inherit2")
    artifact2.save()
    yield artifact1, artifact2
    artifact1.delete(permanent=True, storage=True)
    artifact2.delete(permanent=True, storage=True)


# also see test_feature_manager!
def test_add_from(get_test_artifacts):
    artifact1, artifact2 = get_test_artifacts
    label_names = [f"Project {i}" for i in range(3)]
    records = [ln.Record(name=label_name) for label_name in label_names]
    ln.save(records)

    cell_line_names = [f"Cell line {i}" for i in range(3)]
    cell_lines = [bt.CellLine(name=name) for name in cell_line_names]
    ln.save(cell_lines)

    # pass a list of length 0
    artifact2.labels.add([])
    # now actually pass the labels
    artifact2.labels.add(records)
    # here test add without passing a feature
    artifact2.labels.add(cell_lines)
    assert artifact2.cell_lines.count() == len(cell_lines)

    assert artifact1.records.exists() is False
    artifact1.labels.add_from(artifact2)
    assert artifact1.records.count() == artifact2.records.count()
    assert artifact1.cell_lines.count() == artifact2.cell_lines.count()

    artifact2.cell_lines.remove(*cell_lines)
    artifact1.cell_lines.remove(*cell_lines)
    artifact2.records.remove(*records)
    artifact1.records.remove(*records)

    for record in records:
        record.delete(permanent=True)
    for cell_line in cell_lines:
        cell_line.delete(permanent=True)
