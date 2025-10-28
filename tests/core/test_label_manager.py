from pathlib import Path

import bionty as bt
import lamindb as ln
import pytest
from _dataset_fixtures import (  # noqa
    get_mini_csv,
)
from lamindb.errors import ValidationError
from lamindb.models.artifact import add_labels


@pytest.fixture(scope="module")
def adata():
    adata = ln.examples.datasets.anndata_with_obs()
    # add another column
    adata.obs["cell_type_by_expert"] = adata.obs["cell_type"]
    adata.obs.loc["obs0", "cell_type_by_expert"] = "B cell"
    return adata


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


def test_labels_add_using_anndata(adata):
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


def test_labels_get(get_mini_csv: Path):  # noqa: F811
    artifact = ln.Artifact(get_mini_csv, description="test")
    # feature doesn't exist
    with pytest.raises(TypeError):
        artifact.labels.get("x")  # type: ignore
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
