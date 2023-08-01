import lnschema_bionty as lb
import pytest

import lamindb as ln

adata = ln.dev.datasets.anndata_with_obs()
# add another column
adata.obs["cell_type_from_expert"] = adata.obs["cell_type"]
adata.obs.loc["obs0", "cell_type_from_expert"] = "B cell"


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
    feature = ln.Feature.filter(name="project").one()
    assert feature.type == "category"
    assert feature.registries == "core.Label"
    file.delete(storage=True)


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

    # check the basic construction of the feature set based on obs
    feature_set_obs = file.feature_sets.filter(
        ref_field__startswith="core.Feature", filefeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 5
    assert "species" not in feature_set_obs.features.list("name")

    # now, we add species and run checks
    file.features.add_labels(species, feature="species")
    feature = ln.Feature.filter(name="species").one()
    assert feature.type == "category"
    assert feature.registries == "bionty.Species"
    feature_set_obs = file.feature_sets.filter(
        ref_field__startswith="core.Feature", filefeatureset__slot="obs"
    ).one()
    assert feature_set_obs.n == 5
    feature_set_ext = file.feature_sets.filter(
        ref_field__startswith="core.Feature", filefeatureset__slot="ext"
    ).one()
    assert feature_set_ext.n == 1
    assert "species" in feature_set_ext.features.list("name")

    # now we add cell types & tissues and run checks
    file.features.add_labels(cell_types + cell_types_from_expert)
    file.features.add_labels(tissues, feature="tissue")
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
    file.features.add_labels(diseases, feature="disease")
    df = file.features["obs"].df()
    assert set(df["name"]) == {
        "cell_type",
        "cell_type_id",
        "disease",
        "tissue",
        "cell_type_from_expert",
    }
    assert set(df["type"]) == {
        "category",
    }
    assert set(df["registries"]) == {
        "bionty.CellType",
        None,
        "core.Label",
        "bionty.Tissue|core.Label",
    }

    # now, let's add another feature to ext
    project_1 = ln.Label(name="Project 1")
    file.features.add_labels(project_1, feature="project")
    df = file.features["ext"].df()
    assert set(df["name"]) == {
        "species",
        "project",
    }
    assert set(df["type"]) == {
        "category",
    }

    assert set(file.features.get_labels("project").list("name")) == {"Project 1"}
    assert set(file.features.get_labels("disease").list("name")) == {
        "chronic kidney disease",
        "Alzheimer disease",
        "liver lymphoma",
        "cardiac ventricle disorder",
    }
    assert set(file.features.get_labels("species").list("name")) == {"mouse"}
    assert set(file.features.get_labels("tissue")["bionty.Tissue"].list("name")) == {
        "liver",
        "heart",
        "kidney",
        "brain",
    }
    assert set(file.features.get_labels("tissue")["core.Label"].list("name")) == {
        "organoid",
    }
    # currently, we can't stratify the two cases below
    assert set(file.features.get_labels("cell_type").list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }
    assert set(file.features.get_labels("cell_type_from_expert").list("name")) == {
        "T cell",
        "my new cell type",
        "hepatocyte",
        "hematopoietic stem cell",
        "B cell",
    }

    assert set(df["registries"]) == {"bionty.Species", "core.Label"}
    assert project_1 in file.labels.all()

    # call describe
    file.describe()

    # clean up
    ln.Feature.filter(name="species").one().delete()
    ln.File.filter(description="Mini adata").one().delete(storage=True)
    ln.FeatureSet.filter().all().delete()
