import lnschema_bionty as lb
import pytest

import lamindb as ln

lb.settings.auto_save_parents = False


adata = ln.dev.datasets.anndata_with_obs()
# add another column
adata.obs["cell_type_from_expert"] = adata.obs["cell_type"]
adata.obs.loc["obs0", "cell_type_from_expert"] = "B cell"


def test_labels_add():
    label = ln.ULabel(name="Experiment 1")
    file = ln.File(adata, description="test")
    file.save()
    experiment = ln.Feature(name="experiment", type="category")
    with pytest.raises(ValueError) as error:
        file.labels.add("experiment_1", experiment)
    assert (
        error.exconly()
        == "ValueError: Please pass a record (a `Registry` object), not a string, e.g.,"
        " via: label = ln.ULabel(name='experiment_1')"
    )
    with pytest.raises(ln.dev.exceptions.ValidationError) as error:
        file.labels.add(label, experiment)
    assert "not validated. If it looks correct: record.save()" in error.exconly()
    label.save()
    with pytest.raises(TypeError) as error:
        file.labels.add(label, "experiment 1")
    with pytest.raises(ln.dev.exceptions.ValidationError) as error:
        file.labels.add(label, feature=experiment)
    assert (
        error.exconly()
        == "lamindb.dev.exceptions.ValidationError: Feature not validated. If it looks"
        " correct: ln.Feature(name='experiment', type='category',"
        " registries='core.ULabel').save()"
    )
    experiment.save()
    file.labels.add(label, feature=experiment)
    # check that the feature was updated with registries = "core.ULabel"
    feature = ln.Feature.filter(name="experiment").one()
    assert feature.registries == "core.ULabel"
    with pytest.raises(TypeError):
        experiments = file.labels.get("experiment")
    # check that the label is there, it's exactly one label with name "Experiment 1"
    experiments = file.labels.get(experiment)
    assert experiments.one().name == "Experiment 1"

    # try adding the same label again, nothing should happen
    file.labels.add(label, feature=experiment)
    # check that the label is there, it's exactly one label with name "Experiment 1"
    experiments = file.labels.get(experiment)
    assert experiments.get().name == "Experiment 1"

    feature_set_n1 = ln.FeatureSet.filter(features__name="experiment").one()

    # running from_values to load validated label records under the hood
    experiment = ln.Feature(
        name="experiment_with_reg", type="category", registries=[ln.ULabel]
    )
    experiment.save()
    ln.ULabel(name="Experiment 2").save()
    file.labels.add("Experiment 2", experiment)
    experiments = file.labels.get(experiment)
    assert experiments.get().name == "Experiment 2"

    # now, try adding a new label for a new feature, extending the feature set
    project = ln.ULabel(name="project 1")
    project.save()
    ln.Feature(name="project", type="category").save()
    features = ln.Feature.lookup()
    file.labels.add(project, feature=features.project)
    # check that the label is there, it's exactly one label with name "Experiment 1"
    projects = file.labels.get(features.project)
    assert projects.get().name == "project 1"

    # here, we test that feature_set_n1 was removed because it was no longer
    # linked to any file
    feature_set_n2 = ln.FeatureSet.filter(features__name="experiment").one()
    assert feature_set_n1.id != feature_set_n2.id
    assert file.feature_sets.get() == feature_set_n2

    # test add_from
    dataset = ln.Dataset(file, name="My dataset")
    dataset.save()
    dataset.labels.add_from(file)
    experiments = dataset.labels.get(experiment)
    assert experiments.get().name == "Experiment 2"

    # test features._add_from
    # first, remove all feature sets
    feature_sets = dataset.feature_sets.all()
    for feature_set in feature_sets:
        dataset.feature_sets.remove(feature_set)
    assert len(dataset.feature_sets.all()) == 0
    # second,
    dataset.features._add_from(file)
    assert set(dataset.feature_sets.all()) == set(feature_sets)

    dataset.delete(permanent=True, storage=True)
    ln.Feature.filter().all().delete()
    ln.ULabel.filter().all().delete()
    ln.FeatureSet.filter().all().delete()


def test_add_labels_using_anndata():
    organism = lb.Organism.from_bionty(name="mouse")
    cell_types = [lb.CellType(name=name) for name in adata.obs["cell_type"].unique()]
    ln.save(cell_types)
    inspector = lb.CellType.inspect(adata.obs["cell_type_from_expert"].unique())
    ln.save([lb.CellType(name=name) for name in inspector.non_validated])
    cell_types_from_expert = lb.CellType.from_values(
        adata.obs["cell_type_from_expert"].unique()
    )
    actual_tissues = [lb.Tissue(name=name) for name in adata.obs["tissue"].unique()]
    organoid = ln.ULabel(name="organoid")
    tissues = actual_tissues + [organoid]
    ln.save(tissues)

    lb.settings.organism = "human"
    lb.settings.auto_save_parents = False

    # clean up DB state
    organism_feature = ln.Feature.filter(name="organism").one_or_none()
    if organism_feature is not None:
        organism_feature.delete()
    file = ln.File.filter(description="Mini adata").one_or_none()
    if file is not None:
        file.delete(permanent=True, storage=True)
    ln.FeatureSet.filter().all().delete()

    # try to construct without registering metadata features
    file = ln.File.from_anndata(
        adata, description="Mini adata", field=lb.Gene.ensembl_gene_id
    )
    assert "obs" not in file._feature_sets
    # add feature set without saving file
    feature_name_feature = ln.Feature(
        name="feature name", type="category", registries="core.ULabel"
    )
    feature_name_feature.save()
    feature_set = ln.FeatureSet(features=[feature_name_feature])
    with pytest.raises(ValueError) as error:
        file.features.add_feature_set(feature_set, slot="random")
    assert (
        error.exconly()
        == "ValueError: Please save the file or dataset before adding a feature set!"
    )

    # now register features we want to validate
    # (we are not interested in cell_type_id, here)
    ln.save(
        ln.Feature.from_df(
            adata.obs[["cell_type", "tissue", "cell_type_from_expert", "disease"]]
        )
    )
    file = ln.File.from_anndata(
        adata, description="Mini adata", field=lb.Gene.ensembl_gene_id
    )
    ln.Feature(name="organism", type="category", registries="bionty.Organism").save()
    features = ln.Feature.lookup()
    with pytest.raises(ValueError) as error:
        file.labels.add(organism, feature=features.organism)
    assert (
        error.exconly()
        == "ValueError: Please save the file/dataset before adding a label!"
    )
    file.save()

    # check the basic construction of the feature set based on obs
    feature_set_obs = file.feature_sets.filter(
        registry="core.Feature", filefeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    assert "organism" not in feature_set_obs.features.list("name")

    # now, we add organism and run checks
    with pytest.raises(ln.dev.exceptions.ValidationError):
        file.labels.add(organism, feature=features.organism)
    organism.save()
    file.labels.add(organism, feature=features.organism)
    feature = ln.Feature.filter(name="organism").one()
    assert feature.type == "category"
    assert feature.registries == "bionty.Organism"
    feature_set_obs = file.feature_sets.filter(
        registry="core.Feature", filefeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    feature_set_ext = file.feature_sets.filter(
        registry="core.Feature", filefeatureset__slot="external"
    ).one()
    assert feature_set_ext.n == 1
    assert "organism" in feature_set_ext.features.list("name")

    # now we add cell types & tissues and run checks
    ln.Feature(name="cell_type", type="category").save()
    ln.Feature(name="cell_type_from_expert", type="category").save()
    ln.Feature(name="tissue", type="category").save()
    file.labels.add(cell_types, feature=features.cell_type)
    file.labels.add(cell_types_from_expert, feature=features.cell_type_from_expert)
    file.labels.add(tissues, feature=features.tissue)
    feature = ln.Feature.filter(name="cell_type").one()
    assert feature.registries == "bionty.CellType"
    feature = ln.Feature.filter(name="cell_type_from_expert").one()
    assert feature.registries == "bionty.CellType"
    feature = ln.Feature.filter(name="tissue").one()
    assert feature.registries == "bionty.Tissue|core.ULabel"
    diseases = [ln.ULabel(name=name) for name in adata.obs["disease"].unique()]
    ln.save(diseases)
    file.labels.add(diseases, feature=features.disease)
    df = file.features["obs"].df()
    assert set(df["name"]) == {
        "cell_type",
        "disease",
        "tissue",
        "cell_type_from_expert",
    }
    assert set(df["type"]) == {
        "category",
    }
    assert set(df["registries"]) == {
        "bionty.CellType",
        "core.ULabel",
        "bionty.Tissue|core.ULabel",
    }

    # now, let's add another feature to ext
    experiment_1 = ln.ULabel(name="experiment_1")
    experiment_1.save()
    ln.Feature(name="experiment", type="category").save()
    features = ln.Feature.lookup()
    file.labels.add(experiment_1, feature=features.experiment)
    df = file.features["external"].df()
    assert set(df["name"]) == {
        "organism",
        "experiment",
    }
    assert set(df["type"]) == {
        "category",
    }

    assert set(file.labels.get(features.experiment).list("name")) == {"experiment_1"}
    assert set(file.labels.get(features.disease).list("name")) == {
        "chronic kidney disease",
        "Alzheimer disease",
        "liver lymphoma",
        "cardiac ventricle disorder",
    }
    assert set(file.labels.get(features.organism).list("name")) == {"mouse"}
    assert set(file.labels.get(features.tissue)["bionty.Tissue"].list("name")) == {
        "liver",
        "heart",
        "kidney",
        "brain",
    }
    assert set(file.labels.get(features.tissue)["core.ULabel"].list("name")) == {
        "organoid",
    }
    # currently, we can't stratify the two cases below
    assert set(file.labels.get(features.cell_type).list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(file.labels.get(features.cell_type, flat_names=True)) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(file.labels.get(features.cell_type_from_expert).list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }

    assert set(df["registries"]) == {"bionty.Organism", "core.ULabel"}
    assert experiment_1 in file.ulabels.all()

    # call describe
    file.describe()

    # clean up
    ln.Feature.filter(name="organism").one().delete()
    ln.File.filter(description="Mini adata").one().delete(permanent=True, storage=True)
    ln.FeatureSet.filter().all().delete()
    feature_name_feature.delete()
    lb.CellType.filter().all().delete()
    lb.Tissue.filter().all().delete()
    lb.Disease.filter().all().delete()
    ln.ULabel.filter().all().delete()


def test_labels_get():
    ln.dev.datasets.file_mini_csv()
    file = ln.File("mini.csv", description="test")
    # feature doesn't exist
    with pytest.raises(TypeError):
        file.labels.get("x")
    # no linked labels
    feature_name_feature = ln.Feature(name="feature name", type="category")
    feature_name_feature.save()
    feature_set = ln.FeatureSet(features=[feature_name_feature])
    feature_set.save()
    file.save()
    assert str(file.features) == "no linked features"
    file.features.add_feature_set(feature_set, slot="random")
    assert file.feature_sets["random"] == feature_set
    file.delete(permanent=True, storage=True)
    feature_set.delete()
    feature_name_feature.delete()
