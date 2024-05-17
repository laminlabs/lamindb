import bionty as bt
import lamindb as ln
import pytest
from lamindb.core.exceptions import ValidationError

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
    experiment = ln.Feature(name="experiment", dtype="cat")
    with pytest.raises(ValidationError):
        artifact.features.add({"experiment": "Experiment 1"})
    experiment.save()
    with pytest.raises(ValidationError) as error:
        artifact.features.add({"experiment": "Experiment 1"})
    assert (
        error.exconly()
        == "lamindb.core.exceptions.ValidationError: Label 'Experiment 1' not found in ln.ULabel"
    )
    ln.ULabel(name="Experiment 1").save()
    artifact.features.add({"experiment": "Experiment 1"})
    assert artifact.artifactulabel_set.first().ulabel.name == "Experiment 1"
    temperature = ln.Feature(name="temperature", dtype="cat").save()
    with pytest.raises(TypeError) as error:
        artifact.features.add({"temperature": 27.2})
    assert (
        error.exconly()
        == "TypeError: Value for feature 'temperature' with type 'cat' must be a string or record."
    )
    temperature.dtype = "number"
    temperature.save()
    artifact.features.add({"temperature": 27.2})
    assert artifact.artifactfeaturevalue_set.first().feature_value.value == 27.2

    # delete everything we created
    artifact.delete(permanent=True)
    ln.ULabel.filter().all().delete()
    ln.Feature.filter().all().delete()
    ln.FeatureSet.filter().all().delete()


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
        " correct: ln.Feature(name='experiment', type='cat[core.ULabel]').save()"
    )
    experiment.save()

    # try to pass list of length zero
    artifact.labels.add([], feature=experiment)
    # now pass a single label
    artifact.labels.add(label, feature=experiment)
    # check that the feature was updated with type = "ULabel"
    feature = ln.Feature.filter(name="experiment").one()
    assert feature.dtype == "cat[core.ULabel]"
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

    feature_set_n1 = ln.FeatureSet.filter(features__name="experiment").one()

    # running from_values to load validated label records under the hood
    experiment = ln.Feature(name="experiment_with_reg", dtype="cat[core.ULabel]")
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

    # here, we test that feature_set_n1 was removed because it was no longer
    # linked to any file
    feature_set_n2 = ln.FeatureSet.filter(features__name="experiment").one()
    assert feature_set_n1.id != feature_set_n2.id
    assert artifact.feature_sets.get() == feature_set_n2

    # test add_from
    collection = ln.Collection(artifact, name="My collection")
    collection.save()
    collection.labels.add_from(artifact)
    experiments = collection.labels.get(experiment)
    assert experiments.get().name == "Experiment 2"

    # test features._add_from
    # first, remove all feature sets
    feature_sets = collection.feature_sets.all()
    for feature_set in feature_sets:
        collection.feature_sets.remove(feature_set)
    assert len(collection.feature_sets.all()) == 0
    # second,
    collection.features._add_from(artifact)
    assert set(collection.feature_sets.all()) == set(feature_sets)

    collection.artifacts[0].delete(permanent=True, storage=True)
    collection.delete(permanent=True)
    ln.Feature.filter().all().delete()
    ln.ULabel.filter().all().delete()
    ln.FeatureSet.filter().all().delete()


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
    feature_name_feature = ln.Feature(name="feature name", dtype="cat[core.ULabel]")
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
    artifact.features.add_from_anndata(var_field=bt.Gene.ensembl_gene_id)

    # check the basic construction of the feature set based on obs
    feature_set_obs = artifact.feature_sets.filter(
        registry="Feature", artifactfeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    assert "organism" not in feature_set_obs.features.list("name")

    # now, we add organism and run checks
    features = ln.Feature.lookup()
    with pytest.raises(ln.core.exceptions.ValidationError):
        artifact.labels.add(organism, feature=features.organism)
    organism.save()
    artifact.labels.add(organism, feature=features.organism)
    organism_link = artifact.artifactorganism_set.first()
    assert organism_link.organism.name == "mouse"
    assert organism_link.feature.name == "organism"
    feature = ln.Feature.filter(name="organism").one()
    assert feature.dtype == "cat[bionty.Organism]"
    feature_set_obs = artifact.feature_sets.filter(
        registry="Feature", artifactfeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    feature_set_ext = artifact.feature_sets.filter(
        registry="Feature", artifactfeatureset__slot="external"
    ).one()
    assert feature_set_ext.n == 1
    assert "organism" in feature_set_ext.features.list("name")

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
    assert feature.dtype == "cat[bionty.Tissue|core.ULabel]"
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
        "cat[core.ULabel]",
        "cat[bionty.Tissue|core.ULabel]",
    }

    # now, let's add another feature to ext
    experiment_1 = ln.ULabel(name="experiment_1")
    experiment_1.save()
    ln.Feature(name="experiment", dtype="cat").save()
    features = ln.Feature.lookup()
    artifact.labels.add(experiment_1, feature=features.experiment)
    df = artifact.features["external"].df()
    assert set(df["name"]) == {
        "organism",
        "experiment",
    }
    assert set(df["dtype"]) == {"cat[bionty.Organism]", "cat[core.ULabel]"}

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

    assert set(df["dtype"]) == {"cat[bionty.Organism]", "cat[core.ULabel]"}
    assert experiment_1 in artifact.ulabels.all()

    # call describe
    artifact.describe()

    # clean up
    bt.Gene.filter().all().delete()
    bt.Organism.filter().all().delete()
    ln.Feature.filter(name="organism").one().delete()
    ln.Artifact.filter(description="Mini adata").one().delete(
        permanent=True, storage=True
    )
    ln.FeatureSet.filter().all().delete()
    feature_name_feature.delete()
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
    assert str(artifact.features) == "no linked features"
    artifact.features.add_feature_set(feature_set, slot="random")
    assert artifact.feature_sets["random"] == feature_set
    artifact.delete(permanent=True, storage=True)
    feature_set.delete()
    feature_name_feature.delete()
