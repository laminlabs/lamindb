import shutil

import bionty as bt
import lamindb as ln
import pytest
import tiledbsoma
import tiledbsoma.io
from lamindb.core import datasets


@pytest.fixture(scope="module")
def small_dataset1_schema():
    # Labels
    ln.ULabel.from_values(["DMSO", "IFNG"], create=True).save()
    ln.ULabel.from_values(
        ["Candidate marker study 1", "Candidate marker study 2"], create=True
    ).save()
    bt.CellType.from_values(["B cell", "T cell"], create=True).save()

    # Features
    # observation-level
    cell_medium = ln.Feature(name="cell_medium", dtype="cat[ULabel]").save()
    sample_note = ln.Feature(name="sample_note", dtype="str").save()
    cell_type_by_expert = ln.Feature(
        name="cell_type_by_expert", dtype="cat[bionty.CellType]"
    ).save()
    cell_type_by_model = ln.Feature(
        name="cell_type_by_model", dtype="cat[bionty.CellType]"
    ).save()
    # dataset-level
    ln.Feature(name="temperature", dtype="float").save()
    ln.Feature(name="study", dtype="cat[ULabel]").save()
    ln.Feature(name="date_of_study", dtype="date").save()
    ln.Feature(name="study_note", dtype="str").save()

    # Schema
    schema = ln.Schema(
        name="small_dataset1_obs_level_metadata",
        otype="DataFrame",
        features=[
            cell_medium,
            sample_note,
            cell_type_by_expert,
            cell_type_by_model,
        ],
        coerce_dtype=True,
    ).save()

    yield schema

    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()


@pytest.fixture(scope="module")
def curator_params():
    """Common curator parameters."""
    return {
        "categoricals": {
            "cell_medium": ln.ULabel.name,
            "cell_type_by_expert": bt.CellType.name,
            "cell_type_by_model": bt.CellType.name,
        },
        "organism": "human",
    }


def test_dataframe_curator(small_dataset1_schema):
    """Test DataFrame curator implementation."""

    df, _ = datasets.small_dataset1(otype="DataFrame")
    curator = ln.curators.DataFrameCurator(df, small_dataset1_schema)
    artifact = curator.save_artifact(key="example_datasets/dataset1.parquet")

    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }
    artifact.delete(permanent=True)


def test_anndata_curator(small_dataset1_schema: ln.Schema, curator_params):
    """Test AnnData curator implementation."""

    obs_schema = small_dataset1_schema
    var_schema = ln.Schema(
        name="small_dataset1_var_schema",
        otype="DataFrame",
        itype="bionty.Gene.ensembl_gene_id",
        dtype="num",
    ).save()
    anndata_schema = ln.Schema(
        name="small_dataset1_anndata_schema",
        otype="AnnData",
        components={"obs": obs_schema, "var": var_schema},
    ).save()

    assert anndata_schema.components.get(slot="obs").composite == anndata_schema
    assert anndata_schema.components.get(slot="var").composite == anndata_schema

    describe_output = anndata_schema.describe(return_str=True)
    assert "small_dataset1_anndata_schema" in describe_output
    assert "small_dataset1_obs_level_metadata" in describe_output
    assert "small_dataset1_var_schema" in describe_output

    adata = datasets.small_dataset1(otype="AnnData")
    curator = ln.curators.AnnDataCurator(adata, anndata_schema)
    artifact = curator.save_artifact(key="example_datasets/dataset1.h5ad")
    assert artifact.schema == anndata_schema

    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }

    artifact.delete(permanent=True)
    obs_schema.delete()
    var_schema.delete()
    anndata_schema.delete()


def test_soma_curator(small_dataset1_schema, curator_params):
    """Test SOMA curator implementation."""
    adata = datasets.small_dataset1(otype="AnnData")
    tiledbsoma.io.from_anndata(
        "./small_dataset1.tiledbsoma", adata, measurement_name="RNA"
    )

    curator = ln.Curator.from_tiledbsoma(
        "./small_dataset1.tiledbsoma",
        var_index={"RNA": ("var_id", bt.Gene.ensembl_gene_id)},
        **curator_params,
    )
    artifact = curator.save_artifact(key="example_datasets/dataset1.tiledbsoma")
    artifact.features.add_values(adata.uns)

    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }

    assert artifact._key_is_virtual
    artifact.delete(permanent=True)
    shutil.rmtree("./small_dataset1.tiledbsoma")
