import shutil

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
import tiledbsoma
import tiledbsoma.io
from lamindb.core import datasets
from lamindb.core.types import FieldAttr
from lamindb.errors import InvalidArgument, ValidationError


@pytest.fixture(scope="module")
def small_dataset1_schema():
    # define labels
    perturbation = ln.ULabel(name="Perturbation", is_type=True).save()
    ln.ULabel(name="DMSO", type=perturbation).save()
    ln.ULabel(name="IFNG", type=perturbation).save()
    ln.ULabel(name="ulabel_but_not_perturbation").save()
    ln.ULabel.from_values(["sample1", "sample2", "sample3"], create=True).save()
    bt.CellType.from_source(name="B cell").save()
    bt.CellType.from_source(name="T cell").save()

    # in next iteration for attrs
    ln.Feature(name="temperature", dtype=float).save()
    # ln.Feature(name="experiment", dtype="cat[ULabel]").save()
    # ln.Feature(name="date_of_study", dtype="date").save()
    # ln.Feature(name="study_note", dtype="str").save()

    # define schema
    schema = ln.Schema(
        name="small_dataset1_obs_level_metadata",
        features=[
            ln.Feature(name="perturbation", dtype=perturbation).save(),
            ln.Feature(name="sample_note", dtype=str).save(),
            ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
            ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
        ],
        index=ln.Feature(name="sample_label", dtype=ln.ULabel).save(),
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
        slots={
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


def test_dataframe_curator(small_dataset1_schema: ln.Schema, ccaplog):
    """Test DataFrame curator implementation."""

    # invalid simple dtype (float)
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

    # Wrong subtype
    df = datasets.small_dataset1(otype="DataFrame", with_wrong_subtype=True)
    curator = ln.curators.DataFrameCurator(df, small_dataset1_schema)
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert (
        error.exconly()
        == """lamindb.errors.ValidationError: 1 term not validated in feature 'perturbation': 'ulabel_but_not_perturbation'
    → fix typos, remove non-existent values, or save terms via: curator.cat.add_new_from('perturbation')
    → a valid label for subtype 'Perturbation' has to be one of ['DMSO', 'IFNG']"""
    )

    # Typo
    df = datasets.small_dataset1(otype="DataFrame", with_typo=True)
    curator = ln.curators.DataFrameCurator(df, small_dataset1_schema)
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert (
        error.exconly()
        == """lamindb.errors.ValidationError: 1 term not validated in feature 'perturbation': 'IFNJ'
    → fix typos, remove non-existent values, or save terms via: curator.cat.add_new_from('perturbation')
    → a valid label for subtype 'Perturbation' has to be one of ['DMSO', 'IFNG']"""
    )

    df = datasets.small_dataset1(otype="DataFrame")
    curator = ln.curators.DataFrameCurator(df, small_dataset1_schema)
    artifact = curator.save_artifact(key="examples/dataset1.parquet")

    assert artifact.schema == small_dataset1_schema
    assert artifact.features.slots["columns"].n == 5
    assert set(artifact.features.get_values()["sample_label"]) == {
        "sample1",
        "sample2",
        "sample3",
    }
    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "CD8-positive, alpha-beta T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }

    # a second dataset with missing values
    ln.ULabel.from_values(["sample4", "sample5", "sample6"], create=True).save()
    df = ln.core.datasets.small_dataset2(otype="DataFrame", gene_symbols_in_index=True)
    curator = ln.curators.DataFrameCurator(df, small_dataset1_schema)
    try:
        curator.validate()
    except ln.errors.ValidationError as error:
        assert str(error).startswith("column 'sample_note' not in dataframe")
    curator.standardize()
    curator.validate()

    artifact.delete(permanent=True)


def test_dataframe_curator_index():
    """Test validating a DataFrame index."""
    df = datasets.small_dataset1(otype="DataFrame", with_index_type_mismatch=True)
    feature = ln.Feature(name="test", dtype="str").save()
    schema = ln.Schema(index=feature).save()
    curator = ln.curators.DataFrameCurator(df, schema)
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: expected series 'None' to have type str"
    )

    schema.delete()
    feature.delete()


def test_dataframe_curator_validate_all_annotate_cat(small_dataset1_schema):
    """Do not pass any features."""

    schema = ln.Schema(itype=ln.Feature).save()
    assert schema.flexible
    df = datasets.small_dataset1(otype="DataFrame")
    artifact = ln.Artifact.from_df(
        df, key="examples/dataset1.parquet", schema=schema
    ).save()
    assert set(artifact.features.get_values()["perturbation"]) == {
        "DMSO",
        "IFNG",
    }
    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "CD8-positive, alpha-beta T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }
    artifact.delete(permanent=True)
    schema.delete()


def test_dataframe_curator_validate_all_annotate_cat2(small_dataset1_schema):
    """Combine half-specifying features, half not."""

    schema = ln.Schema(
        itype=ln.Feature,
        features=[ln.Feature.get(name="perturbation")],
        flexible=True,
    ).save()
    assert schema.flexible
    df = datasets.small_dataset1(otype="DataFrame")
    curator = ln.curators.DataFrameCurator(df, schema)
    artifact = curator.save_artifact(key="examples/dataset1.parquet")
    assert set(artifact.features.get_values()["perturbation"]) == {
        "DMSO",
        "IFNG",
    }
    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "CD8-positive, alpha-beta T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }
    artifact.delete(permanent=True)
    schema.delete()


def test_schema_new_genes(ccaplog):
    df = pd.DataFrame(
        index=pd.Index(
            [
                "ENSG00000139618",  # BRCA2
                "ENSG00000141510",  # TP53
                "ENSG00999000001",  # Invalid ID
                "ENSG00999000002",  # Invalid ID
            ],
            name="ensembl",
        )
    )
    feature = ln.Feature(name="ensembl", dtype=bt.Gene.ensembl_gene_id).save()
    schema = ln.Schema(index=feature).save()
    curator = ln.curators.DataFrameCurator(df, schema)
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: 2 terms not validated in feature 'index': 'ENSG00999000001', 'ENSG00999000002'"
    )

    assert (
        "2 terms not validated in feature 'index': 'ENSG00999000001', 'ENSG00999000002'"
        in ccaplog.text
    )

    schema.delete()
    feature.delete()


def test_schema_no_match_ensembl():
    df = pd.DataFrame(
        index=pd.Index(
            [
                "ENSG99999999998",  # Invalid ID
                "ENSG99999999999",  # Invalid ID
            ],
            name="ensembl",
        )
    )
    schema = ln.Schema(
        index=ln.Feature(name="ensembl", dtype=bt.Gene.ensembl_gene_id).save()
    ).save()
    curator = ln.curators.DataFrameCurator(df, schema)
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert (
        error.exconly()
        == """lamindb.errors.ValidationError: 2 terms not validated in feature 'index': 'ENSG99999999998', 'ENSG99999999999'
    → fix typos, remove non-existent values, or save terms via: curator.cat.add_new_from('index')"""
    )

    schema.delete()


def test_schema_mixed_ensembl_symbols(ccaplog):
    """Quite some datasets have mixed ensembl gene IDs and symbols.

    The expected behavior is that an error is raised when such a dataset is encountered because
    currently lamindb does not support validating values against a union of Fields.

    The current behavior is that these cases automatically pass.
    """
    df = pd.DataFrame(
        index=pd.Index(
            [
                "ENSG00000139618",
                "ENSG00000141510",
                "BRCA2",  # symbol
                "TP53",  # symbol
            ],
            name="ensembl",
        )
    )
    schema = ln.Schema(
        index=ln.Feature(name="ensembl", dtype=bt.Gene.ensembl_gene_id).save()
    ).save()
    curator = ln.curators.DataFrameCurator(df, schema)
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert error.exconly().startswith(
        "lamindb.errors.ValidationError: 2 terms not validated in feature 'index': 'BRCA2', 'TP53'"
    )

    assert "2 terms not validated in feature 'index': 'BRCA2', 'TP53'" in ccaplog.text

    schema.delete()


def test_anndata_curator_different_components(small_dataset1_schema: ln.Schema):
    obs_schema = small_dataset1_schema

    for add_comp in ["var.T", "obs", "uns"]:
        var_schema = ln.Schema(
            name="scRNA_seq_var_schema",
            itype=bt.Gene.ensembl_gene_id,
            dtype="num",
        ).save()

        # always assume var
        components = {"var.T": var_schema}
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
            slots=components,
        ).save()
        assert small_dataset1_schema.id is not None, small_dataset1_schema
        assert anndata_schema.slots["var.T"] == var_schema
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
        assert isinstance(curator.slots["var.T"], ln.curators.DataFrameCurator)
        if add_comp == "obs":
            assert isinstance(curator.slots["obs"], ln.curators.DataFrameCurator)
        if add_comp == "uns":
            assert isinstance(curator.slots["uns"], ln.curators.DataFrameCurator)
        artifact = ln.Artifact.from_anndata(
            adata, key="examples/dataset1.h5ad", schema=anndata_schema
        )
        assert artifact._curator._is_validated  # important test, do not remove
        artifact.save()
        assert not hasattr(artifact, "_curator")  # test that curator is deleted
        assert artifact.schema == anndata_schema
        assert artifact.features.slots["var.T"].n == 3  # 3 genes get linked
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


def test_anndata_curator_varT_curation():
    ln.Schema.filter(itype="bionty.Gene.ensembl_gene_id").delete()
    varT_schema = ln.Schema(itype=bt.Gene.ensembl_gene_id, maximal_set=True).save()
    slot = "var.T"
    components = {slot: varT_schema}
    anndata_schema = ln.Schema(
        otype="AnnData",
        slots=components,
    ).save()
    for with_gene_typo in [True, False]:
        adata = datasets.small_dataset1(otype="AnnData", with_gene_typo=with_gene_typo)
        if with_gene_typo:
            with pytest.raises(ValidationError) as error:
                artifact = ln.Artifact.from_anndata(
                    adata, key="examples/dataset1.h5ad", schema=anndata_schema
                ).save()
            assert error.exconly() == (
                f"lamindb.errors.ValidationError: 1 term not validated in feature 'columns' in slot '{slot}': 'GeneTypo'\n"
                f"    → fix typos, remove non-existent values, or save terms via: curator.slots['{slot}'].cat.add_new_from('columns')"
            )
        else:
            for n_max_records in [2, 4]:
                ln.settings.annotation.n_max_records = n_max_records
                artifact = ln.Artifact.from_anndata(
                    adata, key="examples/dataset1.h5ad", schema=anndata_schema
                ).save()
                assert artifact.features.slots[slot].n == 3  # 3 genes get linked
                assert (
                    artifact.features.slots[slot].itype == "bionty.Gene.ensembl_gene_id"
                )
                if n_max_records == 2:
                    assert not artifact.features.slots[slot].members.exists()
                else:
                    assert artifact.features.slots[slot].members.df()[
                        "ensembl_gene_id"
                    ].tolist() == [
                        "ENSG00000153563",
                        "ENSG00000010610",
                        "ENSG00000170458",
                    ]

                artifact.delete(permanent=True)

            anndata_schema.delete()
            varT_schema.delete()


def test_anndata_curator_varT_curation_legacy(ccaplog):
    varT_schema = ln.Schema(itype=bt.Gene.ensembl_gene_id, maximal_set=True).save()
    slot = "var"
    components = {slot: varT_schema}
    anndata_schema = ln.Schema(
        otype="AnnData",
        slots=components,
    ).save()
    for with_gene_typo in [True, False]:
        adata = datasets.small_dataset1(otype="AnnData", with_gene_typo=with_gene_typo)
        if with_gene_typo:
            with pytest.raises(ValidationError) as error:
                artifact = ln.Artifact.from_anndata(
                    adata, key="examples/dataset1.h5ad", schema=anndata_schema
                ).save()
            assert error.exconly() == (
                f"lamindb.errors.ValidationError: 1 term not validated in feature 'var_index' in slot '{slot}': 'GeneTypo'\n"
                f"    → fix typos, remove non-existent values, or save terms via: curator.slots['{slot}'].cat.add_new_from('var_index')"
            )
        else:
            artifact = ln.Artifact.from_anndata(
                adata, key="examples/dataset1.h5ad", schema=anndata_schema
            ).save()
            assert (
                "auto-transposed `var` for backward compat, please indicate transposition in the schema definition by calling out `.T`: slots={'var.T': itype=bt.Gene.ensembl_gene_id}"
                in ccaplog.text
            )
            assert artifact.features.slots[slot].n == 3  # 3 genes get linked
            assert set(
                artifact.features.slots[slot].members.df()["ensembl_gene_id"]
            ) == {
                "ENSG00000153563",
                "ENSG00000010610",
                "ENSG00000170458",
            }

            artifact.delete(permanent=True)

            anndata_schema.delete()
            varT_schema.delete()


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
    artifact = curator.save_artifact(key="examples/dataset1.tiledbsoma")

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
        slots={"obs": small_dataset1_schema},
    ).save()
    assert small_dataset1_schema.id is not None, small_dataset1_schema
    adata = datasets.small_dataset1(otype="AnnData")
    curator = ln.curators.AnnDataCurator(adata, anndata_schema_no_var)

    artifact = curator.save_artifact(key="examples/dataset1_no_var.h5ad")
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
    ln.settings.verbosity = "hint"
    curator.validate()
    curator.slots["rna:var"].cat.standardize("columns")
    curator.slots["rna:var"].cat.add_new_from("columns")
    artifact = curator.save_artifact(key="mudata_papalexi21_subset.h5mu")
    assert artifact.schema == mudata_schema
    assert set(artifact.features.slots.keys()) == {
        "obs",
        "rna:var",
        "rna:obs",
        "hto:obs",
    }

    artifact.delete(permanent=True)


@pytest.fixture(scope="module")
def spatialdata_blobs_schema():
    import sys
    from pathlib import Path

    docs_path = Path.cwd() / "docs" / "scripts"
    sys.path.append(str(docs_path))

    from define_schema_spatialdata import (
        attrs_schema,
        obs_schema,
        sample_schema,
        tech_schema,
        varT_schema,
    )

    spatialdata_schema_legacy = ln.Schema(
        otype="SpatialData",
        slots={
            "bio": sample_schema,
            "table:obs": obs_schema,
            "table:var": varT_schema,
        },
    ).save()

    spatialdata_schema_new = ln.Schema(
        otype="SpatialData",
        slots={
            "attrs:sample": sample_schema,
            "attrs:tech": tech_schema,
            "attrs": attrs_schema,
            "table:obs": obs_schema,
            "table:var.T": varT_schema,
        },
    ).save()

    yield spatialdata_schema_legacy, spatialdata_schema_new

    from lamindb.models import SchemaComponent

    SchemaComponent.filter().delete()
    spatialdata_schema_legacy.delete()
    spatialdata_schema_new.delete()
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter(type__isnull=False).delete()
    ln.ULabel.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.DevelopmentalStage.filter().delete()
    bt.Disease.filter().delete()


def test_spatialdata_curator(
    spatialdata_blobs_schema: ln.Schema,
):
    spatialdata_schema_legacy, spatialdata_schema_new = spatialdata_blobs_schema
    spatialdata = ln.core.datasets.spatialdata_blobs()

    # wrong dataset
    with pytest.raises(InvalidArgument):
        ln.curators.SpatialDataCurator(pd.DataFrame(), spatialdata_blobs_schema)
    # wrong schema
    with pytest.raises(InvalidArgument):
        ln.curators.SpatialDataCurator(
            spatialdata, spatialdata_schema_legacy.slots["bio"]
        )

    curator = ln.curators.SpatialDataCurator(spatialdata, spatialdata_schema_legacy)
    with pytest.raises(ln.errors.ValidationError):
        curator.validate()
    spatialdata.tables["table"].var.drop(index="ENSG00000999999", inplace=True)

    artifact = ln.Artifact.from_spatialdata(
        spatialdata,
        key="examples/spatialdata1.zarr",
        schema=spatialdata_schema_legacy,
    ).save()
    assert artifact.schema == spatialdata_schema_legacy
    assert artifact.features.slots.keys() == {
        "bio",
        "table:var",
        "table:obs",
    }
    assert artifact.features.get_values()["disease"] == "Alzheimer disease"
    artifact.delete(permanent=True)

    artifact = ln.Artifact.from_spatialdata(
        spatialdata,
        key="examples/spatialdata1.zarr",
        schema=spatialdata_schema_new,
    ).save()
    assert artifact.schema == spatialdata_schema_new
    assert artifact.features.slots.keys() == {
        "attrs:bio",
        "attrs:tech",
        "attrs",
        "tables:table:obs",
        "tables:table:var.T",
    }
    assert artifact.features.get_values()["assay"] == "Visium Spatial Gene Expression"
    assert (
        artifact.features.describe(return_str=True)
        == """Artifact .zarr/SpatialData
└── Dataset features
    ├── attrs:bio • 2       [Feature]
    │   developmental_sta…  cat[bionty.Devel…  adult stage
    │   disease             cat[bionty.Disea…  Alzheimer disease
    ├── attrs:tech • 1      [Feature]
    │   assay               cat[bionty.Exper…  Visium Spatial Gene Expression
    ├── attrs • 2           [Feature]
    │   bio                 dict
    │   tech                dict
    ├── tables:table:obs …  [Feature]
    │   sample_region       str
    └── tables:table:var.…  [bionty.Gene.ens…
        BRCA2               num
        BRAF                num"""
    )

    artifact.delete(permanent=True)


def test_tiledbsoma_curator(small_dataset1_schema: ln.Schema, clean_soma_files):
    """Test TiledbTiledbsomaExperimentCurator with schema."""
    obs_schema = ln.Schema(
        features=[
            ln.Feature(name="cell_type_by_expert", dtype=bt.CellType).save(),
            ln.Feature(name="cell_type_by_model", dtype=bt.CellType).save(),
        ],
    ).save()

    var_schema = ln.Schema(
        features=[
            ln.Feature(name="var_id", dtype=bt.Gene.ensembl_gene_id).save(),
        ],
        coerce_dtype=True,
    ).save()

    soma_schema = ln.Schema(
        otype="tiledbsoma",
        slots={
            "obs": obs_schema,
            "ms:RNA": var_schema,
        },
    ).save()

    # Convert AnnData to SOMA format
    adata = ln.core.datasets.small_dataset1(otype="AnnData")
    tiledbsoma.io.from_anndata(
        "small_dataset.tiledbsoma", adata, measurement_name="RNA"
    )
    experiment = tiledbsoma.Experiment.open("small_dataset.tiledbsoma")

    # Test with invalid dataset
    with pytest.raises(ln.errors.InvalidArgument):
        ln.curators.TiledbsomaExperimentCurator(adata, soma_schema)

    # Test with invalid schema
    with pytest.raises(ln.errors.InvalidArgument):
        ln.curators.TiledbsomaExperimentCurator(experiment, small_dataset1_schema)

    curator = ln.curators.TiledbsomaExperimentCurator(experiment, soma_schema)

    assert "obs" in curator.slots
    assert "ms:RNA" in curator.slots

    curator.validate()

    artifact = curator.save_artifact(
        key="examples/soma_experiment.tiledbsoma",
        description="SOMA experiment with schema validation",
    )

    assert artifact.schema == soma_schema
    assert "obs" in artifact.features.slots
    assert "ms:RNA" in artifact.features.slots

    # Check feature values are properly annotated
    assert set(artifact.features.get_values()["cell_type_by_expert"]) == {
        "CD8-positive, alpha-beta T cell",
        "B cell",
    }
    assert set(artifact.features.get_values()["cell_type_by_model"]) == {
        "T cell",
        "B cell",
    }

    # Altered data (gene typo)
    adata_typo = ln.core.datasets.small_dataset1(otype="AnnData", with_gene_typo=True)
    typo_soma_path = "./small_dataset1_typo.tiledbsoma"
    tiledbsoma.io.from_anndata(typo_soma_path, adata_typo, measurement_name="RNA")
    experiment_typo = tiledbsoma.Experiment.open(typo_soma_path)

    curator_typo = ln.curators.TiledbsomaExperimentCurator(experiment_typo, soma_schema)

    # Validation should fail due to typo
    with pytest.raises(ln.errors.ValidationError) as error:
        curator_typo.validate()
    assert "GeneTypo" in str(error.value)

    # Clean up
    artifact.delete(permanent=True)
    soma_schema.delete()
    var_schema.delete()
    obs_schema.delete()
