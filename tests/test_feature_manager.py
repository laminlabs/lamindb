import lnschema_bionty as lb
import pytest

import lamindb as ln

lb.settings.auto_save_parents = False


adata = ln.dev.datasets.anndata_with_obs()
# add another column
adata.obs["cell_type_from_expert"] = adata.obs["cell_type"]
adata.obs.loc["obs0", "cell_type_from_expert"] = "B cell"


def test_add_labels():
    label = ln.Label(name="Experiment 1")
    file = ln.File(adata, description="test")
    file.save()
    experiment = ln.Feature(name="experiment", type="category")
    with pytest.raises(ValueError) as error:
        file.add_labels("experiment_1", experiment)
    assert (
        error.exconly()
        == "ValueError: Please pass a record (a `Registry` object), not a string, e.g.,"
        " via: label = ln.Label(name='experiment_1')"
    )
    with pytest.raises(ln.dev.exceptions.ValidationError) as error:
        file.add_labels(label, experiment)
    assert "not validated. If it looks correct: record.save()" in error.exconly()
    label.save()
    with pytest.raises(TypeError) as error:
        file.add_labels(label, "experiment 1")
    with pytest.raises(ln.dev.exceptions.ValidationError) as error:
        file.add_labels(label, feature=experiment)
    assert (
        error.exconly()
        == "lamindb.dev.exceptions.ValidationError: Feature not validated. If it looks"
        " correct: ln.Feature(name='experiment', type='category',"
        " registries='core.Label').save()"
    )
    experiment.save()
    file.add_labels(label, feature=experiment)
    # check that the feature was updated with registries = "core.Label"
    feature = ln.Feature.filter(name="experiment").one()
    assert feature.registries == "core.Label"
    with pytest.raises(TypeError):
        experiments = file.get_labels("experiment")
    # check that the label is there, it's exactly one label with name "Experiment 1"
    experiments = file.get_labels(experiment)
    assert experiments.get().name == "Experiment 1"

    # try adding the same label again, nothing should happen
    file.add_labels(label, feature=experiment)
    # check that the label is there, it's exactly one label with name "Experiment 1"
    experiments = file.get_labels(experiment)
    assert experiments.get().name == "Experiment 1"

    feature_set_n1 = ln.FeatureSet.filter(features__name="experiment").one()

    # now, try adding a new label for a new feature, extending the feature set
    project = ln.Label(name="project 1")
    project.save()
    ln.Feature(name="project", type="category").save()
    features = ln.Feature.lookup()
    file.add_labels(project, feature=features.project)
    # check that the label is there, it's exactly one label with name "Experiment 1"
    projects = file.get_labels(features.project)
    assert projects.get().name == "project 1"

    # here, we test that feature_set_n1 was removed because it was no longer
    # linked to any file
    feature_set_n2 = ln.FeatureSet.filter(features__name="experiment").one()
    assert feature_set_n1.id != feature_set_n2.id
    assert file.feature_sets.get() == feature_set_n2

    file.delete(storage=True)
    ln.Feature.filter().all().delete()
    ln.Label.filter().all().delete()
    ln.FeatureSet.filter().all().delete()


def test_add_labels_using_anndata():
    species = lb.Species.from_bionty(name="mouse")
    cell_types = [lb.CellType(name=name) for name in adata.obs["cell_type"].unique()]
    ln.save(cell_types)
    cell_types_from_expert = [
        lb.CellType(name=name) for name in adata.obs["cell_type_from_expert"].unique()
    ]
    ln.save(cell_types_from_expert)
    actual_tissues = [lb.Tissue(name=name) for name in adata.obs["tissue"].unique()]
    organoid = ln.Label(name="organoid")
    tissues = actual_tissues + [organoid]
    ln.save(tissues)

    lb.settings.species = "human"
    lb.settings.auto_save_parents = False

    # clean up DB state
    species_feature = ln.Feature.filter(name="species").one_or_none()
    if species_feature is not None:
        species_feature.delete()
    file = ln.File.filter(description="Mini adata").one_or_none()
    if file is not None:
        file.delete(storage=True)
    ln.FeatureSet.filter().all().delete()

    # try to construct without registering metadata features
    file = ln.File.from_anndata(
        adata, description="Mini adata", field=lb.Gene.ensembl_gene_id
    )
    assert "obs" not in file._feature_sets
    # add feature set without saving file
    feature_name_feature = ln.Feature(
        name="feature name", type="category", registries="core.Label"
    )
    feature_name_feature.save()
    feature_set = ln.FeatureSet(features=[feature_name_feature])
    with pytest.raises(ValueError):
        file.features.add_feature_set(feature_set, slot="random")

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
    ln.Feature(name="species", type="category", registries="bionty.Species").save()
    features = ln.Feature.lookup()
    with pytest.raises(ValueError) as error:
        file.add_labels(species, feature=features.species)
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
    assert "species" not in feature_set_obs.features.list("name")

    # now, we add species and run checks
    with pytest.raises(ln.dev.exceptions.ValidationError):
        file.add_labels(species, feature=features.species)
    species.save()
    file.add_labels(species, feature=features.species)
    feature = ln.Feature.filter(name="species").one()
    assert feature.type == "category"
    assert feature.registries == "bionty.Species"
    feature_set_obs = file.feature_sets.filter(
        registry="core.Feature", filefeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    feature_set_ext = file.feature_sets.filter(
        registry="core.Feature", filefeatureset__slot="external"
    ).one()
    assert feature_set_ext.n == 1
    assert "species" in feature_set_ext.features.list("name")

    # now we add cell types & tissues and run checks
    ln.Feature(name="cell_type", type="category").save()
    ln.Feature(name="cell_type_from_expert", type="category").save()
    ln.Feature(name="tissue", type="category").save()
    file.add_labels(cell_types, feature=features.cell_type)
    file.add_labels(cell_types_from_expert, feature=features.cell_type_from_expert)
    file.add_labels(tissues, feature=features.tissue)
    feature = ln.Feature.filter(name="cell_type").one()
    assert feature.registries == "bionty.CellType"
    feature = ln.Feature.filter(name="cell_type_from_expert").one()
    assert feature.registries == "bionty.CellType"
    feature = ln.Feature.filter(name="tissue").one()
    assert feature.registries == "bionty.Tissue|core.Label"
    diseases = [ln.Label(name=name) for name in adata.obs["disease"].unique()]
    ln.save(diseases)
    file.add_labels(diseases, feature=features.disease)
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
        "core.Label",
        "bionty.Tissue|core.Label",
    }

    # now, let's add another feature to ext
    experiment_1 = ln.Label(name="experiment_1")
    experiment_1.save()
    ln.Feature(name="experiment", type="category").save()
    features = ln.Feature.lookup()
    file.add_labels(experiment_1, feature=features.experiment)
    df = file.features["external"].df()
    assert set(df["name"]) == {
        "species",
        "experiment",
    }
    assert set(df["type"]) == {
        "category",
    }

    assert set(file.get_labels(features.experiment).list("name")) == {"experiment_1"}
    assert set(file.get_labels(features.disease).list("name")) == {
        "chronic kidney disease",
        "Alzheimer disease",
        "liver lymphoma",
        "cardiac ventricle disorder",
    }
    assert set(file.get_labels(features.species).list("name")) == {"mouse"}
    assert set(file.get_labels(features.tissue)["bionty.Tissue"].list("name")) == {
        "liver",
        "heart",
        "kidney",
        "brain",
    }
    assert set(file.get_labels(features.tissue)["core.Label"].list("name")) == {
        "organoid",
    }
    # currently, we can't stratify the two cases below
    assert set(file.get_labels(features.cell_type).list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(file.get_labels(features.cell_type, flat_names=True)) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(file.get_labels(features.cell_type_from_expert).list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }

    assert set(df["registries"]) == {"bionty.Species", "core.Label"}
    assert experiment_1 in file.labels.all()

    # call describe
    file.describe()

    # clean up
    ln.Feature.filter(name="species").one().delete()
    ln.File.filter(description="Mini adata").one().delete(storage=True)
    ln.FeatureSet.filter().all().delete()
    feature_name_feature.delete()
    lb.CellType.filter().all().delete()
    lb.Tissue.filter().all().delete()
    lb.Disease.filter().all().delete()
    ln.Label.filter().all().delete()


def test_get_labels():
    ln.dev.datasets.file_mini_csv()
    file = ln.File("mini.csv", description="test")
    # feature doesn't exist
    with pytest.raises(TypeError):
        file.get_labels("x")
    # no linked labels
    feature_name_feature = ln.Feature(name="feature name", type="category")
    feature_name_feature.save()
    feature_set = ln.FeatureSet(features=[feature_name_feature])
    feature_set.save()
    file.save()
    assert str(file.features) == "no linked features"
    file.features.add_feature_set(feature_set, slot="random")
    assert file.features.get_feature_set(slot="random") == feature_set
    file.delete(storage=True)
    feature_set.delete()
    feature_name_feature.delete()
