import bionty as bt
import lamindb as ln
import pytest
from lamindb.core.exceptions import ValidationError
from lnschema_core.models import FeatureValue

bt.settings.auto_save_parents = False


@pytest.fixture(scope="module")
def adata():
    adata = ln.core.datasets.anndata_with_obs()
    # add another column
    adata.obs["cell_type_from_expert"] = adata.obs["cell_type"]
    adata.obs.loc["obs0", "cell_type_from_expert"] = "B cell"
    return adata


# below is the main new test for the main way of annotating with
# features
def test_features_add(adata):
    ln.ULabel(name="Experiment 1")
    artifact = ln.Artifact.from_anndata(adata, description="test")
    artifact.save()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"experiment": "Experiment 1"})
    assert error.exconly().startswith(
        "lamindb.core.exceptions.ValidationError: These keys could not be validated:"
    )
    experiment = ln.Feature(name="experiment", dtype="cat")
    experiment.save()
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values({"experiment": "Experiment 1"})
    assert error.exconly().startswith(
        "lamindb.core.exceptions.ValidationError: These values could not be validated: ['Experiment 1']"
    )
    experiment_label = ln.ULabel(name="Experiment 1").save()
    # add the label without the feature first
    artifact.ulabels.add(experiment_label)
    assert artifact.ulabel_links.get().ulabel.name == "Experiment 1"
    assert artifact.ulabel_links.get().feature is None

    # now add the label with the feature and make sure that it has the feature annotation
    artifact.features.add_values({"experiment": "Experiment 1"})
    assert artifact.ulabel_links.get().ulabel.name == "Experiment 1"
    assert artifact.ulabel_links.get().feature.name == "experiment"
    # repeat
    artifact.features.add_values({"experiment": "Experiment 1"})
    assert artifact.ulabel_links.get().ulabel.name == "Experiment 1"

    # numerical feature
    temperature = ln.Feature(name="temperature", dtype="cat").save()
    with pytest.raises(TypeError) as error:
        artifact.features.add_values({"temperature": 27.2})
    assert (
        error.exconly()
        == "TypeError: Value for feature 'temperature' with type 'cat' must be a string or record."
    )
    temperature.dtype = "number"
    temperature.save()
    artifact.features.add_values({"temperature": 27.2})
    assert artifact.feature_values.first().value == 27.2

    features = {
        "experiment": "Experiment 2",
        "project": "project_1",
        "is_validated": True,
        "cell_type_by_expert": "T Cell",
        "temperature": 100.0,
        "donor": "U0123",
    }
    with pytest.raises(ValidationError) as error:
        artifact.features.add_values(features)
    print(error.exconly())
    assert (
        error.exconly()
        == """\
lamindb.core.exceptions.ValidationError: These keys could not be validated: ['project', 'is_validated', 'cell_type_by_expert', 'donor']
If there are no typos, create features for them:

  ln.Feature(name='project', dtype='cat[ULabel]').save()
  ln.Feature(name='is_validated', dtype='bool').save()
  ln.Feature(name='cell_type_by_expert', dtype='cat[ULabel]').save()
  ln.Feature(name='donor', dtype='cat[ULabel]').save()"""
    )

    ln.Feature(name="project", dtype="cat[ULabel]").save()
    ln.Feature(name="is_validated", dtype="bool").save()
    ln.Feature(name="cell_type_by_expert", dtype="cat[ULabel]").save()
    ln.Feature(name="donor", dtype="cat[ULabel]").save()

    with pytest.raises(ValidationError) as error:
        artifact.features.add_values(features)
    print(error.exconly())
    assert (
        error.exconly()
        == """\
lamindb.core.exceptions.ValidationError: These values could not be validated: ['Experiment 2', 'project_1', 'T Cell', 'U0123']
If there are no typos, create ulabels for them:

  ulabels = ln.ULabel.from_values(['Experiment 2', 'project_1', 'T Cell', 'U0123'], create=True)
  ln.save(ulabels)"""
    )

    ulabels = ln.ULabel.from_values(
        ["Experiment 2", "project_1", "T Cell", "U0123"], create=True
    )
    ln.save(ulabels)

    artifact.features.add_values(features)
    assert set(artifact.feature_values.all().values_list("value", flat=True)) == {
        27.2,
        True,
        100.0,
    }

    assert ln.Artifact.filter(feature_values__value=27.2).one()

    print(artifact.features.get_values())
    print(artifact.features.__repr__())
    # hard to test because of italic formatting
    msg = """\
    'experiment' = 'Experiment 1', 'Experiment 2'
    'project' = 'project_1'
    'cell_type_by_expert' = 'T Cell'
    'donor' = 'U0123'
    'is_validated' = True
    'temperature' = 27.2, 100.0
"""
    assert artifact.features.__repr__().endswith(msg)
    assert artifact.features.get_values() == {
        "experiment": ["Experiment 1", "Experiment 2"],
        "project": "project_1",
        "cell_type_by_expert": "T Cell",
        "donor": "U0123",
        "is_validated": True,
        "temperature": [27.2, 100.0],
    }

    # repeat
    artifact.features.add_values(features)
    assert set(artifact.feature_values.all().values_list("value", flat=True)) == {
        27.2,
        True,
        100.0,
    }
    assert artifact.features.__repr__().endswith(msg)

    with pytest.raises(ValidationError) as error:
        ln.Artifact.features.filter(
            temperature_with_typo=100.0, project="project_1"
        ).one()
    assert error.exconly().startswith(
        "lamindb.core.exceptions.ValidationError: Some keys in the filter expression are not registered as features:"
    )

    ln.Artifact.features.filter(temperature=100.0).one()
    ln.Artifact.features.filter(project="project_1").one()
    ln.Artifact.features.filter(is_validated=True).one()
    ln.Artifact.features.filter(
        temperature=100.0, project="project_1", donor="U0123"
    ).one()

    # delete everything we created
    artifact.delete(permanent=True)
    ln.ULabel.filter().all().delete()
    ln.FeatureSet.filter().all().delete()
    ln.Feature.filter().all().delete()


def test_labels_add(adata):
    label = ln.ULabel(name="Experiment 1")
    artifact = ln.Artifact.from_anndata(adata, description="test")
    artifact.save()
    experiment = ln.Feature(name="experiment", dtype="cat")
    with pytest.raises(ValueError) as error:
        artifact.labels.add("experiment_1", experiment)
    assert (
        error.exconly()
        == "ValueError: Please pass a record (a `Registry` object), not a string, e.g.,"
        " via: label = ln.ULabel(name='experiment_1')"
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
        == "lamindb.core.exceptions.ValidationError: Feature not validated. If it looks"
        " correct: ln.Feature(name='experiment', type='cat[ULabel]').save()"
    )
    experiment.save()

    # try to pass list of length zero
    artifact.labels.add([], feature=experiment)
    # now pass a single label
    artifact.labels.add(label, feature=experiment)
    # check that the feature was updated with type = "ULabel"
    feature = ln.Feature.filter(name="experiment").one()
    assert feature.dtype == "cat[ULabel]"
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
    experiment = ln.Feature(name="experiment_with_reg", dtype="cat[ULabel]")
    experiment.save()
    ln.ULabel(name="Experiment 2").save()
    artifact.labels.add("Experiment 2", experiment)
    experiments = artifact.labels.get(experiment)
    assert experiments.get().name == "Experiment 2"

    # now, try adding a new label for a new feature, extending the feature set
    project = ln.ULabel(name="project 1")
    project.save()
    ln.Feature(name="project", dtype="cat").save()
    features = ln.Feature.lookup()
    artifact.labels.add(project, feature=features.project)
    # check that the label is there, it's exactly one label with name "Experiment 1"
    projects = artifact.labels.get(features.project)
    assert projects.get().name == "project 1"

    # test add_from
    collection = ln.Collection(artifact, name="My collection")
    collection.save()
    from lamindb.core._label_manager import get_labels_as_dict

    collection.labels.add_from(artifact)
    experiments = collection.labels.get(experiment)
    assert experiments.get().name == "Experiment 2"

    collection.delete(permanent=True)
    artifact.delete(permanent=True)
    ln.FeatureSet.filter().all().delete()
    ln.Feature.filter().all().delete()
    ln.ULabel.filter().all().delete()


def test_add_labels_using_anndata(adata):
    organism = bt.Organism.from_public(name="mouse")
    cell_types = [bt.CellType(name=name) for name in adata.obs["cell_type"].unique()]
    ln.save(cell_types)
    inspector = bt.CellType.inspect(adata.obs["cell_type_from_expert"].unique())
    ln.save([bt.CellType(name=name) for name in inspector.non_validated])
    cell_types_from_expert = bt.CellType.from_values(
        adata.obs["cell_type_from_expert"].unique()
    )
    actual_tissues = [bt.Tissue(name=name) for name in adata.obs["tissue"].unique()]
    organoid = ln.ULabel(name="organoid")
    tissues = actual_tissues + [organoid]
    ln.save(tissues)

    bt.settings.organism = "human"
    bt.settings.auto_save_parents = False

    # clean up DB state
    organism_feature = ln.Feature.filter(name="organism").one_or_none()
    if organism_feature is not None:
        organism_feature.delete()
    artifact = ln.Artifact.filter(description="Mini adata").one_or_none()
    if artifact is not None:
        artifact.delete(permanent=True, storage=True)
    ln.FeatureSet.filter().all().delete()

    # try to construct without registering metadata features
    artifact = ln.Artifact.from_anndata(adata, description="Mini adata")
    if not artifact._state.adding:
        artifact.delete(permanent=True)  # make sure we get a fresh one
        artifact = ln.Artifact.from_anndata(adata, description="Mini adata")
    # add feature set without saving file
    feature_name_feature = ln.Feature(name="feature name", dtype="cat[ULabel]")
    feature_name_feature.save()
    feature_set = ln.FeatureSet(features=[feature_name_feature])
    with pytest.raises(ValueError) as error:
        artifact.features.add_feature_set(feature_set, slot="random")
    assert (
        error.exconly()
        == "ValueError: Please save the artifact or collection before adding a feature"
        " set!"
    )

    # now register features we want to validate
    # (we are not interested in cell_type_id, here)
    ln.save(
        ln.Feature.from_df(
            adata.obs[["cell_type", "tissue", "cell_type_from_expert", "disease"]]
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
    feature_set_obs = artifact.feature_sets.filter(
        registry="Feature", artifact_links__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    assert "organism" not in feature_set_obs.features.list("name")

    # now, we add organism and run checks
    features = ln.Feature.lookup()
    with pytest.raises(ln.core.exceptions.ValidationError):
        artifact.labels.add(organism, feature=features.organism)
    organism.save()
    artifact.labels.add(organism, feature=features.organism)
    organism_link = artifact.organism_links.first()
    assert organism_link.organism.name == "mouse"
    assert organism_link.feature.name == "organism"
    feature = ln.Feature.filter(name="organism").one()
    assert feature.dtype == "cat[bionty.Organism]"
    feature_set_obs = artifact.feature_sets.filter(
        registry="Feature", artifact_links__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    # TODO, write a test that queries the organism feature
    # assert "organism" in feature_set_ext.features.list("name")

    # now we add cell types & tissues and run checks
    ln.Feature(name="cell_type", dtype="cat").save()
    ln.Feature(name="cell_type_from_expert", dtype="cat").save()
    ln.Feature(name="tissue", dtype="cat").save()
    artifact.labels.add(cell_types, feature=features.cell_type)
    artifact.labels.add(cell_types_from_expert, feature=features.cell_type_from_expert)
    artifact.labels.add(tissues, feature=features.tissue)
    feature = ln.Feature.filter(name="cell_type").one()
    assert feature.dtype == "cat[bionty.CellType]"
    feature = ln.Feature.filter(name="cell_type_from_expert").one()
    assert feature.dtype == "cat[bionty.CellType]"
    feature = ln.Feature.filter(name="tissue").one()
    assert feature.dtype == "cat[bionty.Tissue|ULabel]"
    diseases = [ln.ULabel(name=name) for name in adata.obs["disease"].unique()]
    ln.save(diseases)
    artifact.labels.add(diseases, feature=features.disease)
    df = artifact.features["obs"].df()
    assert set(df["name"]) == {
        "cell_type",
        "disease",
        "tissue",
        "cell_type_from_expert",
    }
    assert set(df["dtype"]) == {
        "cat[bionty.CellType]",
        "cat[ULabel]",
        "cat[bionty.Tissue|ULabel]",
    }

    # now, let's add another feature to ext
    experiment_1 = ln.ULabel(name="experiment_1")
    experiment_1.save()
    ln.Feature(name="experiment", dtype="cat").save()
    features = ln.Feature.lookup()
    artifact.labels.add(experiment_1, feature=features.experiment)
    # TODO: replace the following with an updated test
    # df = artifact.features["external"].df()
    # assert set(df["name"]) == {
    #     "organism",
    #     "experiment",
    # }
    # assert set(df["dtype"]) == {"cat[bionty.Organism]", "cat[ULabel]"}

    assert set(artifact.labels.get(features.experiment).list("name")) == {
        "experiment_1"
    }
    assert set(artifact.labels.get(features.disease).list("name")) == {
        "chronic kidney disease",
        "Alzheimer disease",
        "liver lymphoma",
        "cardiac ventricle disorder",
    }
    assert set(artifact.labels.get(features.organism).list("name")) == {"mouse"}
    assert set(artifact.labels.get(features.tissue)["bionty.Tissue"].list("name")) == {
        "liver",
        "heart",
        "kidney",
        "brain",
    }
    assert set(artifact.labels.get(features.tissue)["ULabel"].list("name")) == {
        "organoid",
    }
    # currently, we can't stratify the two cases below
    assert set(artifact.labels.get(features.cell_type).list("name")) == {
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
    assert set(artifact.labels.get(features.cell_type_from_expert).list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert experiment_1 in artifact.ulabels.all()

    # call describe
    artifact.describe()

    # clean up
    artifact.delete(permanent=True)
    bt.Gene.filter().all().delete()
    bt.Organism.filter().all().delete()
    ln.FeatureSet.filter().all().delete()
    ln.Feature.filter().all().delete()
    bt.CellType.filter().all().delete()
    bt.Tissue.filter().all().delete()
    bt.Disease.filter().all().delete()
    ln.ULabel.filter().all().delete()


def test_labels_get():
    ln.core.datasets.file_mini_csv()
    artifact = ln.Artifact("mini.csv", description="test")
    # feature doesn't exist
    with pytest.raises(TypeError):
        artifact.labels.get("x")
    # no linked labels
    feature_name_feature = ln.Feature(name="feature name", dtype="cat")
    feature_name_feature.save()
    feature_set = ln.FeatureSet(features=[feature_name_feature])
    feature_set.save()
    artifact.save()
    assert str(artifact.features) == ""
    artifact.features.add_feature_set(feature_set, slot="random")
    assert artifact.feature_sets.first() == feature_set
    artifact.delete(permanent=True, storage=True)
    feature_set.delete()
    feature_name_feature.delete()
