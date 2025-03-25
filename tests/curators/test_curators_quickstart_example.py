import shutil

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
import tiledbsoma
import tiledbsoma.io
from lamindb.core import datasets
from lamindb.core.types import FieldAttr
from lamindb.errors import InvalidArgument


@pytest.fixture(scope="module")
def small_dataset1_schema():
    # define labels
    perturbation = ln.ULabel(name="Perturbation", is_type=True).save()
    ln.ULabel(name="DMSO", type=perturbation).save()
    ln.ULabel(name="IFNG", type=perturbation).save()
    bt.CellType.from_source(name="B cell").save()
    bt.CellType.from_source(name="T cell").save()

    # in next iteration for attrs
    ln.Feature(name="temperature", dtype=float).save()
    # ln.Feature(name="study", dtype="cat[ULabel]").save()
    # ln.Feature(name="date_of_study", dtype="date").save()
    # ln.Feature(name="study_note", dtype="str").save()

    # define schema
    schema = ln.Schema(
        name="small_dataset1_obs_level_metadata",
        features=[
            ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
            ln.Feature(name="sample_note", dtype=str).save(),
            ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
            ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
        ],
    ).save()

    yield schema

    schema.delete()
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter(type__isnull=False).delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()


@pytest.fixture(scope="module")
def curator_params():
    """Common curator parameters."""
    return {
        "categoricals": {
            "perturbation": ln.ULabel.name,
            "cell_type_by_expert": bt.CellType.name,
            "cell_type_by_model": bt.CellType.name,
        },
        "organism": "human",
    }


@pytest.fixture(scope="module")
def mudata_papalexi21_subset_schema():
    # define labels
    perturbation = ln.ULabel(name="Perturbation", is_type=True).save()
    ln.ULabel(name="Perturbed", type=perturbation).save()
    ln.ULabel(name="NT", type=perturbation).save()

    replicate = ln.ULabel(name="Replicate", is_type=True).save()
    ln.ULabel(name="rep1", type=replicate).save()
    ln.ULabel(name="rep2", type=replicate).save()
    ln.ULabel(name="rep3", type=replicate).save()

    # define obs schema
    obs_schema = ln.Schema(
        name="mudata_papalexi21_subset_obs_schema",
        features=[
            ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
            ln.Feature(name="replicate", dtype="cat[ULabel[Replicate]]").save(),
        ],
    ).save()

    obs_schema_rna = ln.Schema(
        name="mudata_papalexi21_subset_rna_obs_schema",
        features=[
            ln.Feature(name="nCount_RNA", dtype=int).save(),
            ln.Feature(name="nFeature_RNA", dtype=int).save(),
            ln.Feature(name="percent.mito", dtype=float).save(),
        ],
        coerce_dtype=True,
    ).save()

    obs_schema_hto = ln.Schema(
        name="mudata_papalexi21_subset_hto_obs_schema",
        features=[
            ln.Feature(name="nCount_HTO", dtype=float).save(),
            ln.Feature(name="nFeature_HTO", dtype=int).save(),
            ln.Feature(name="technique", dtype=bt.ExperimentalFactor).save(),
        ],
        coerce_dtype=True,
    ).save()

    var_schema_rna = ln.Schema(
        name="mudata_papalexi21_subset_rna_var_schema",
        itype=bt.Gene.symbol,
        dtype=float,
    ).save()

    # define composite schema
    mudata_schema = ln.Schema(
        name="mudata_papalexi21_subset_mudata_schema",
        otype="MuData",
        components={
            "obs": obs_schema,
            "rna:obs": obs_schema_rna,
            "hto:obs": obs_schema_hto,
            "rna:var": var_schema_rna,
        },
    ).save()

    yield mudata_schema

    mudata_schema.delete()
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter(type__isnull=False).delete()
    ln.ULabel.filter().delete()
    bt.ExperimentalFactor.filter().delete()


@pytest.fixture(scope="module")
def spatialdata_blobs_schema():
    sample_schema = ln.Schema(
        name="blobs_sample_level_metadata",
        features=[
            ln.Feature(name="assay", dtype=bt.ExperimentalFactor).save(),
            ln.Feature(name="disease", dtype=bt.Disease).save(),
            ln.Feature(name="developmental_stage", dtype=bt.DevelopmentalStage).save(),
        ],
        coerce_dtype=True,
    ).save()

    blobs_obs_schema = ln.Schema(
        name="blobs_obs_level_metadata",
        features=[
            ln.Feature(name="sample_region", dtype="str").save(),
        ],
        coerce_dtype=True,
    ).save()

    blobs_var_schema = ln.Schema(
        name="visium_var_schema", itype=bt.Gene.ensembl_gene_id, dtype=int
    ).save()

    spatialdata_schema = ln.Schema(
        name="blobs_spatialdata_schema",
        otype="SpatialData",
        components={
            "sample": sample_schema,
            "table:obs": blobs_obs_schema,
            "table:var": blobs_var_schema,
        },
    ).save()

    yield spatialdata_schema

    from lamindb.models import SchemaComponent

    SchemaComponent.filter().delete()
    spatialdata_schema.delete()
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter(type__isnull=False).delete()
    ln.ULabel.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.DevelopmentalStage.filter().delete()
    bt.Disease.filter().delete()


def test_dataframe_curator(small_dataset1_schema: ln.Schema):
    """Test DataFrame curator implementation."""

    feature_to_fail = ln.Feature(name="treatment_time_h", dtype=float).save()
    schema = ln.Schema(
        name="small_dataset1_obs_level_metadata_v2",
        features=[
            ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]").save(),
            ln.Feature(name="sample_note", dtype=str).save(),
            ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
            ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
            feature_to_fail,
        ],
    ).save()

    df = datasets.small_dataset1(otype="DataFrame")
    curator = ln.curators.DataFrameCurator(df, schema)
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert (
        error.exconly()
        == "lamindb.errors.ValidationError: Column 'treatment_time_h' failed series or dataframe validator 0: <Check check_function: Column 'treatment_time_h' failed dtype check for 'float': got int64>"
    )
    schema.delete()
    feature_to_fail.delete()
    curator = ln.curators.DataFrameCurator(df, small_dataset1_schema)
    artifact = curator.save_artifact(key="example_datasets/dataset1.parquet")

    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "CD8-positive, alpha-beta T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }

    # a second dataset with missing values
    df = ln.core.datasets.small_dataset2(otype="DataFrame", gene_symbols_in_index=True)
    curator = ln.curators.DataFrameCurator(df, small_dataset1_schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert str(error).startswith("column 'sample_note' not in dataframe")
    curator.standardize()
    curator.validate()

    artifact.delete(permanent=True)


def test_anndata_curator(small_dataset1_schema: ln.Schema):
    """Test AnnData curator implementation."""

    obs_schema = small_dataset1_schema

    for add_comp in ["var", "obs", "uns"]:
        var_schema = ln.Schema(
            name="scRNA_seq_var_schema",
            itype=bt.Gene.ensembl_gene_id,
            dtype="num",
        ).save()

        # always assume var
        components = {"var": var_schema}
        if add_comp == "obs":
            components["obs"] = obs_schema
        if add_comp == "uns":
            uns_schema = ln.Schema(
                name="flexible_uns_schema",
                itype=ln.Feature,
            ).save()
            components["uns"] = uns_schema

        anndata_schema = ln.Schema(
            name="small_dataset1_anndata_schema",
            otype="AnnData",
            components=components,
        ).save()
        assert small_dataset1_schema.id is not None, small_dataset1_schema
        assert anndata_schema.slots["var"] == var_schema
        if add_comp == "obs":
            assert anndata_schema.slots["obs"] == obs_schema
        if add_comp == "uns":
            assert anndata_schema.slots["uns"] == uns_schema

        describe_output = anndata_schema.describe(return_str=True)
        assert "small_dataset1_anndata_schema" in describe_output
        assert "scRNA_seq_var_schema" in describe_output
        if add_comp == "obs":
            assert "small_dataset1_obs_level_metadata" in describe_output
        if add_comp == "uns":
            assert "flexible_uns_schema" in describe_output

        adata = datasets.small_dataset1(otype="AnnData")
        curator = ln.curators.AnnDataCurator(adata, anndata_schema)
        assert isinstance(curator.slots["var"], ln.curators.DataFrameCurator)
        if add_comp == "obs":
            assert isinstance(curator.slots["obs"], ln.curators.DataFrameCurator)
        if add_comp == "uns":
            assert isinstance(curator.slots["uns"], ln.curators.DataFrameCurator)
        artifact = curator.save_artifact(key="example_datasets/dataset1.h5ad")
        assert artifact.schema == anndata_schema
        assert artifact.features.slots["var"].n == 3  # 3 genes get linked
        if add_comp == "obs":
            assert artifact.features.slots["obs"] == obs_schema
            # deprecated
            assert artifact.features._schema_by_slot["obs"] == obs_schema
            assert artifact.features._feature_set_by_slot["obs"] == obs_schema

            assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
                "CD8-positive, alpha-beta T cell",
                "B cell",
            }
            assert set(artifact.features.get_values()["cell_type_by_model"]) == {
                "T cell",
                "B cell",
            }
        if add_comp == "uns":
            assert artifact.features.slots["uns"].features.first() == ln.Feature.get(
                name="temperature"
            )

        artifact.delete(permanent=True)
        anndata_schema.delete()
        var_schema.delete()


def test_soma_curator(
    small_dataset1_schema: ln.Schema, curator_params: dict[str, str | FieldAttr]
):
    """Test SOMA curator implementation."""
    adata = datasets.small_dataset1(otype="AnnData")
    tiledbsoma.io.from_anndata(
        "./small_dataset1.tiledbsoma", adata, measurement_name="RNA"
    )

    curator = ln.Curator.from_tiledbsoma(  # type: ignore
        "./small_dataset1.tiledbsoma",
        var_index={"RNA": ("var_id", bt.Gene.ensembl_gene_id)},
        **curator_params,
    )
    artifact = curator.save_artifact(key="example_datasets/dataset1.tiledbsoma")
    # artifact.features.add_values(adata.uns)

    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "CD8-positive, alpha-beta T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }

    assert artifact._key_is_virtual
    artifact.delete(permanent=True)
    shutil.rmtree("./small_dataset1.tiledbsoma")


def test_anndata_curator_no_var(small_dataset1_schema: ln.Schema):
    assert small_dataset1_schema.id is not None, small_dataset1_schema
    # test no var schema
    anndata_schema_no_var = ln.Schema(
        name="small_dataset1_anndata_schema_no_var",
        otype="AnnData",
        components={"obs": small_dataset1_schema},
    ).save()
    assert small_dataset1_schema.id is not None, small_dataset1_schema
    adata = datasets.small_dataset1(otype="AnnData")
    curator = ln.curators.AnnDataCurator(adata, anndata_schema_no_var)
    artifact = curator.save_artifact(key="example_datasets/dataset1_no_var.h5ad")
    artifact.delete(permanent=True)
    anndata_schema_no_var.delete()


def test_mudata_curator(
    mudata_papalexi21_subset_schema: ln.Schema, small_dataset1_schema: ln.Schema
):
    mudata_schema = mudata_papalexi21_subset_schema
    mdata = ln.core.datasets.mudata_papalexi21_subset()
    # TODO: refactor organism
    bt.settings.organism = "human"
    # wrong dataset
    with pytest.raises(InvalidArgument):
        ln.curators.MuDataCurator(pd.DataFrame(), mudata_schema)
    # wrong schema
    with pytest.raises(InvalidArgument):
        ln.curators.MuDataCurator(mdata, small_dataset1_schema)
    curator = ln.curators.MuDataCurator(mdata, mudata_schema)
    assert curator.slots.keys() == {
        "obs",
        "rna:obs",
        "hto:obs",
        "rna:var",
    }
    artifact = curator.save_artifact(key="mudata_papalexi21_subset.h5mu")
    assert artifact.schema == mudata_schema
    assert artifact.features.slots.keys() == {
        "obs",
        "['rna'].var",
        "['rna'].obs",
        "['hto'].obs",
    }

    artifact.delete(permanent=True)


def test_spatialdata_curator(
    spatialdata_blobs_schema: ln.Schema, small_dataset1_schema: ln.Schema
):
    spatialdata_schema = spatialdata_blobs_schema
    spatialdata = ln.core.datasets.spatialdata_blobs()

    # wrong dataset
    with pytest.raises(InvalidArgument):
        ln.curators.SpatialDataCurator(pd.DataFrame(), spatialdata_blobs_schema)
    # wrong schema
    with pytest.raises(InvalidArgument):
        ln.curators.SpatialDataCurator(spatialdata, small_dataset1_schema)

    curator = ln.curators.SpatialDataCurator(spatialdata, spatialdata_schema)
    try:
        curator.validate()
    except ln.errors.ValidationError:
        pass

    # validate again (must pass now) and save artifact
    artifact = curator.save_artifact(key="example_datasets/spatialdata1.zarr")
    assert artifact.schema == spatialdata_schema
    assert artifact.features.slots.keys() == {
        "sample",
        "['table'].var",
        "['table'].obs",
    }

    artifact.delete(permanent=True)
