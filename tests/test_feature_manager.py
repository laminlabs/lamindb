import lnschema_bionty as lb
import pytest

import lamindb as ln
from lamindb.dev._feature_manager import create_features_df

lb.settings.auto_save_parents = False


adata = ln.dev.datasets.anndata_with_obs()
# add another column
adata.obs["cell_type_from_expert"] = adata.obs["cell_type"]
adata.obs.loc["obs0", "cell_type_from_expert"] = "B cell"


def test_features_add_labels():
    label = ln.Label(name="Experiment 1")
    file = ln.File(adata, description="test")
    file.save()
    with pytest.raises(ln.dev.exc.ValidationError) as error:
        file.add_labels(label)
    assert "not validated. If it looks correct: record.save()" in error.exconly()
    label.save()
    with pytest.raises(ValueError) as error:
        file.add_labels(label)
    assert (
        error.exconly()
        == "ValueError: Please pass feature: add_labels(labels, feature='myfeature')"
    )
    with pytest.raises(ln.dev.exc.ValidationError) as error:
        file.add_labels(label, feature="experiment")
    assert (
        error.exconly()
        == "lamindb.dev.exc.ValidationError: Feature not validated. If it looks"
        " correct: ln.Feature(name='experiment', type='category',"
        " registries='core.Label').save()"
    )
    ln.Feature(name="experiment", type="category", registries="core.Label").save()
    file.add_labels(label, feature="experiment")
    file.add_labels(ln.Label.filter(name="Experiment 1").all(), feature="experiment")
    feature = ln.Feature.filter(name="experiment").one()
    assert feature.type == "category"
    assert feature.registries == "core.Label"

    file.delete(storage=True)
    ln.Feature.filter().all().delete()
    ln.Label.filter().all().delete()
    ln.FeatureSet.filter().all().delete()


def test_features_add_labels_using_anndata():
    species = lb.Species.from_bionty(name="mouse")
    cell_types = lb.CellType.from_values(adata.obs["cell_type"], "name")
    ln.save(cell_types)
    cell_types_from_expert = lb.CellType.from_values(
        adata.obs["cell_type_from_expert"], "name"
    )
    ln.save(cell_types_from_expert)
    actual_tissues = lb.Tissue.from_values(adata.obs["tissue"], "name")
    organoid = ln.Label(name="organoid")
    tissues = actual_tissues + [organoid]
    ln.save(tissues)

    assert cell_types[0]._feature == "cell_type"
    assert cell_types[-1]._feature == "cell_type"
    assert tissues[0]._feature == "tissue"

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
        adata, description="Mini adata", var_ref=lb.Gene.ensembl_gene_id
    )
    assert "obs" not in file._feature_sets
    # add feature set without saving file
    feature = ln.Feature(name="feature name", type="category", registries="core.Label")
    feature_set = ln.FeatureSet(features=[feature])
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
        adata, description="Mini adata", var_ref=lb.Gene.ensembl_gene_id
    )

    with pytest.raises(ValueError) as error:
        file.add_labels("species")
    assert (
        error.exconly()
        == "ValueError: Please pass a record (a `Registry` object), not a string, e.g.,"
        " via: label = ln.Label(name='species')"
    )

    with pytest.raises(ValueError) as error:
        file.add_labels(species, feature="species")
    assert (
        error.exconly()
        == "ValueError: Please save the file or dataset before adding a label!"
    )

    file.save()

    # check the basic construction of the feature set based on obs
    feature_set_obs = file.feature_sets.filter(
        registry="core.Feature", filefeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    assert "species" not in feature_set_obs.features.list("name")

    # now, we add species and run checks
    with pytest.raises(ln.dev.exc.ValidationError):
        file.add_labels(species, feature="species")
    species.save()
    with pytest.raises(ln.dev.exc.ValidationError):
        file.add_labels(species, feature="species")
    ln.Feature(name="species", type="category", registries="bionty.Species").save()
    file.add_labels(species, feature="species")
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
    file.add_labels(cell_types + cell_types_from_expert)
    file.add_labels(tissues, feature="tissue")
    feature = ln.Feature.filter(name="cell_type").one()
    assert feature.type == "category"
    assert feature.registries == "bionty.CellType"
    feature = ln.Feature.filter(name="cell_type_from_expert").one()
    assert feature.type == "category"
    assert feature.registries == "bionty.CellType"
    feature = ln.Feature.filter(name="tissue").one()
    assert feature.type == "category"
    assert feature.registries == "bionty.Tissue|core.Label"
    diseases = ln.Label.from_values(adata.obs["disease"])
    ln.save(diseases)
    file.add_labels(diseases, feature="disease")
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
    file.add_labels(experiment_1, feature="experiment")
    df = file.features["external"].df()
    assert set(df["name"]) == {
        "species",
        "experiment",
    }
    assert set(df["type"]) == {
        "category",
    }

    assert set(file.get_labels("experiment").list("name")) == {"experiment_1"}
    assert set(file.get_labels("disease").list("name")) == {
        "chronic kidney disease",
        "Alzheimer disease",
        "liver lymphoma",
        "cardiac ventricle disorder",
    }
    assert set(file.get_labels("species").list("name")) == {"mouse"}
    assert set(file.get_labels("tissue")["bionty.Tissue"].list("name")) == {
        "liver",
        "heart",
        "kidney",
        "brain",
    }
    assert set(file.get_labels("tissue")["core.Label"].list("name")) == {
        "organoid",
    }
    # currently, we can't stratify the two cases below
    assert set(file.get_labels("cell_type").list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(file.get_labels("cell_type", flat_names=True)) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(file.get_labels("cell_type_from_expert").list("name")) == {
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


def test_get_labels():
    ln.dev.datasets.file_mini_csv()
    file = ln.File("mini.csv", description="test")
    # feature doesn't exist
    with pytest.raises(ValueError):
        file.get_labels("x")

    # no linked labels
    feature = ln.Feature(name="feature name", type="category")
    feature_set = ln.FeatureSet(features=[feature])
    feature_set.save()
    file.save()
    assert str(file.features) == "no linked features"
    file.features.add_feature_set(feature_set, slot="random")
    with pytest.raises(ValueError):
        file.get_labels("feature name")

    # exclude=False for create_features_df
    df = create_features_df(file, [feature_set], exclude=False)
    assert df.shape[0] == 1
    file.delete(storage=True)
    feature_set.delete()
