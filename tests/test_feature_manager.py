import lnschema_bionty as lb
import pytest

import lamindb as ln

adata = ln.dev.datasets.anndata_with_obs()


def test_features_add_labels():
    label = ln.Label(name="Project 1")
    label.save()
    file = ln.File(adata)
    file.save()
    with pytest.raises(ValueError) as error:
        file.features.add_labels(label)
    assert (
        error.exconly()
        == "ValueError: Please pass feature: add_labels(labels, feature='myfeature')"
    )
    file.features.add_labels(label, feature="project")
    feature = ln.Feature.select(name="project").one()
    assert feature.type == "category"
    assert feature.label_orms == "core.Label"
    file.delete(storage=True)


def test_features_add_labels_using_anndata():
    species = lb.Species.from_bionty(name="mouse")
    cell_types = lb.CellType.from_values(adata.obs["cell_type"], "name")
    tissues = lb.Tissue.from_values(adata.obs["tissue"], "name")

    assert species._feature == "species"
    assert cell_types[0]._feature == "cell_type"
    assert cell_types[-1]._feature == "cell_type"
    assert tissues[0]._feature == "tissue"

    lb.settings.species = "human"
    lb.settings.auto_save_parents = False

    # clean up DB state
    species_feature = ln.Feature.select(name="species").one_or_none()
    if species_feature is not None:
        species_feature.delete()
    file = ln.File.select(description="Mini adata").one_or_none()
    if file is not None:
        file.delete(storage=True)
    ln.FeatureSet.select().all().delete()

    file = ln.File.from_anndata(
        adata, description="Mini adata", var_ref=lb.Gene.ensembl_gene_id
    )

    with pytest.raises(ValueError) as error:
        file.features.add_labels("species")
    assert (
        error.exconly()
        == "ValueError: Please pass a record (an ORM object), not a string, e.g., via:"
        " label = ln.Label(name='species')"
    )

    with pytest.raises(ValueError) as error:
        file.features.add_labels(species, feature="species")
    assert (
        error.exconly()
        == "ValueError: Please save the file or dataset before adding a label!"
    )

    file.save()

    feature_set_obs = file.feature_sets.filter(
        ref_field__startswith="core.Feature", filefeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    assert "species" not in feature_set_obs.features.list("name")

    file.features.add_labels(species, feature="species")
    feature = ln.Feature.select(name="species").one()
    assert feature.type == "category"
    assert feature.label_orms == "bionty.Species"

    feature_set_obs = file.feature_sets.filter(
        ref_field__startswith="core.Feature", filefeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 4
    feature_set_ext = file.feature_sets.filter(
        ref_field__startswith="core.Feature", filefeatureset__slot="ext"
    ).one()
    assert feature_set_ext.n == 1
    assert "species" in feature_set_ext.features.list("name")

    file.features.add_labels(cell_types + tissues)
    feature = ln.Feature.select(name="cell_type").one()
    assert feature.type == "category"
    assert feature.label_orms == "bionty.CellType"
    feature = ln.Feature.select(name="tissue").one()
    assert feature.type == "category"
    assert feature.label_orms == "bionty.Tissue"

    # on purpose, we don't use bionty ORM here, to simulate an ordinary label
    diseases = ln.Label.from_values(adata.obs["disease"])
    file.features.add_labels(diseases, feature="disease")

    df = file.features["obs"].df()
    assert set(df["name"]) == {
        "cell_type",
        "cell_type_id",
        "disease",
        "tissue",
    }
    assert set(df["type"]) == {
        "category",
    }
    assert set(df["label_orms"]) == {
        "bionty.CellType",
        None,
        "core.Label",
        "bionty.Tissue",
    }

    df = file.features["ext"].df()
    assert set(df["name"]) == {
        "species",
    }
    assert set(df["type"]) == {
        "category",
    }
    assert set(df["label_orms"]) == {"bionty.Species"}

    # clean up
    ln.Feature.select(name="species").one().delete()
    ln.File.select(description="Mini adata").one().delete(storage=True)
    ln.FeatureSet.select().all().delete()
