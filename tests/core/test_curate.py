from unittest.mock import Mock

import anndata as ad
import bionty as bt
import lamindb as ln
import mudata as md
import pandas as pd
import pytest
from lamindb._curate import CurateLookup, ValidationError


@pytest.fixture(scope="module")
def df():
    return pd.DataFrame(
        {
            "cell_type": [
                "cerebral pyramidal neuron",  # on purpose, should be "cerebral cortex pyramidal neuron"
                "astrocytic glia",  # synonym of astrocyte
                "oligodendrocyte",
            ],
            "cell_type_2": ["oligodendrocyte", "oligodendrocyte", "astrocyte"],
            "assay_ontology_id": ["EFO:0008913", "EFO:0008913", "EFO:0008913"],
            "donor": ["D0001", "D0002", "DOOO3"],
        }
    )


@pytest.fixture(scope="module")
def adata():
    df = pd.DataFrame(
        {
            "cell_type": [
                "cerebral cortex pyramidal neuron",
                "astrocytic glia",  # synonym of astrocyte
                "oligodendrocyte",
            ],
            "cell_type_2": [
                "oligodendrocyte",
                "oligodendrocyte",
                "astrocyte",
            ],
            "assay_ontology_id": ["EFO:0008913", "EFO:0008913", "EFO:0008913"],
            "donor": ["D0001", "D0002", "DOOO3"],
        }
    )
    df.index = ["obs1", "obs2", "obs3"]

    X = pd.DataFrame(
        {
            "TCF-1": [1, 2, 3],  # synonym of TCF7
            "PDCD1": [4, 5, 6],
            "CD3E": [7, 8, 9],
            "CD4": [10, 11, 12],
            "CD8A": [13, 14, 15],
        },
        index=["obs1", "obs2", "obs3"],
    )

    return ad.AnnData(X=X, obs=df)


@pytest.fixture(scope="module")
def mdata(adata):
    mdata = md.MuData({"rna": adata, "rna_2": adata})
    mdata.obs["donor"] = ["D0001", "D0002", "DOOO3"]

    return mdata


@pytest.fixture(scope="module")
def categoricals():
    return {
        "cell_type": bt.CellType.name,
        "cell_type_2": bt.CellType.name,
        "assay_ontology_id": bt.ExperimentalFactor.ontology_id,
        "donor": ln.ULabel.name,
    }


@pytest.fixture
def curate_lookup(categoricals):
    return CurateLookup(categoricals=categoricals, using_key="undefined")


@pytest.fixture
def mock_registry():
    registry = Mock()
    registry.lookup = Mock(return_value="mocked lookup")
    return registry


@pytest.fixture
def mock_transform():
    mock_transform = ln.Transform(name="mock", version="0.0.0", type="notebook")
    mock_transform.save()
    return mock_transform


def test_df_curator(df, categoricals):
    df = df.copy()
    curator = ln.Curator.from_df(df, categoricals=categoricals)
    with pytest.raises(ValidationError):
        _ = curator.non_validated
    validated = curator.validate()
    assert curator.non_validated == {
        "cell_type": ["cerebral pyramidal neuron", "astrocytic glia"],
        "donor": ["D0001", "D0002", "DOOO3"],
    }
    assert validated is False

    # deprecated method
    curator.add_new_from_columns()

    # standardize
    with pytest.raises(KeyError):
        curator.standardize("nonexistent-key")
    curator.standardize("all")
    assert curator.non_validated == {
        "cell_type": ["cerebral pyramidal neuron"],
        "donor": ["D0001", "D0002", "DOOO3"],
    }
    assert "astrocyte" in df["cell_type"].values

    # add new
    curator.add_new_from("donor")
    assert curator.non_validated == {"cell_type": ["cerebral pyramidal neuron"]}

    # lookup
    cell_types = curator.lookup(public=True)["cell_type"]
    df["cell_type"] = df["cell_type"].replace(
        {"cerebral pyramidal neuron": cell_types.cerebral_cortex_pyramidal_neuron.name}
    )
    validated = curator.validate()
    assert validated is True
    assert curator.non_validated == {}

    # no need to standardize
    curator.standardize("cell_type")

    artifact = curator.save_artifact(description="test-curate-df")

    assert (
        artifact.cell_types.through.filter(artifact_id=artifact.id)
        .df()["label_ref_is_name"]
        .values.sum()
        == 5
    )
    assert (
        artifact.cell_types.through.filter(artifact_id=artifact.id)
        .df()["feature_ref_is_name"]
        .values.sum()
        == 5
    )

    assert (
        artifact.experimental_factors.through.filter(artifact_id=artifact.id)
        .df()["label_ref_is_name"]
        .values.sum()
        == 0
    )
    assert (
        artifact.experimental_factors.through.filter(artifact_id=artifact.id)
        .df()["feature_ref_is_name"]
        .values.sum()
        == 1
    )

    assert set(artifact.features.get_values()["cell_type"]) == {
        "cerebral cortex pyramidal neuron",
        "astrocyte",
        "oligodendrocyte",
    }
    assert set(artifact.features.get_values()["cell_type_2"]) == {
        "oligodendrocyte",
        "astrocyte",
    }

    # clean up
    artifact.delete(permanent=True)
    ln.ULabel.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.CellType.filter().delete()
    ln.FeatureSet.filter().delete()


def test_custom_using_invalid_field_lookup(curate_lookup):
    with pytest.raises(AttributeError) as excinfo:
        _ = curate_lookup["invalid_field"]
    assert '"CurateLookup" object has no attribute "invalid_field"' in str(
        excinfo.value
    )


def test_additional_args_with_all_key(df, categoricals):
    curator = ln.Curator.from_df(df, categoricals=categoricals)
    with pytest.raises(ValueError) as error:
        curator.add_new_from("all", extra_arg="not_allowed")
    assert "Cannot pass additional arguments to 'all' key!" in str(error.value)


def test_save_columns_not_defined_in_fields(df, categoricals):
    curator = ln.Curator.from_df(df, categoricals=categoricals)
    with pytest.raises(ValidationError) as error:
        curator._update_registry("nonexistent")
    assert "Feature nonexistent is not part of the fields!" in str(error.value)


def test_unvalidated_data_object(df, categoricals):
    curator = ln.Curator.from_df(df, categoricals=categoricals)
    with pytest.raises(ValidationError) as error:
        curator.save_artifact()
    assert "Dataset does not validate. Please curate." in str(error.value)


def test_clean_up_failed_runs():
    mock_transform = ln.Transform()
    mock_transform.save()
    mock_run = ln.Run(mock_transform)
    mock_run.save()
    mock_run_2 = ln.Run(mock_transform)
    mock_run_2.save()

    # Set the default currently used transform and mock run -> these should not be cleaned up
    from lamindb.core._context import context

    previous_transform = context._transform
    previous_run = context.run

    context._transform = mock_transform
    context._run = mock_run

    assert len(ln.Run.filter(transform=mock_transform).all()) == 2

    curator = ln.Curator.from_df(pd.DataFrame())
    curator.clean_up_failed_runs()

    assert len(ln.Run.filter(transform=mock_transform).all()) == 1

    # Revert to old run context to not infer with tests that need the run context
    context._transform = previous_transform
    context._run = previous_run


@pytest.mark.parametrize("to_add", ["donor", "all"])
def test_anndata_curator(adata, categoricals, to_add):
    adata = adata.copy()
    # # must pass an organism
    # with pytest.raises(ValidationError):
    #     ln.Curator.from_anndata(
    #         adata,
    #         categoricals=categoricals,
    #         var_index=bt.Gene.symbol,
    #     ).validate()

    curator = ln.Curator.from_anndata(
        adata,
        categoricals=categoricals,
        var_index=bt.Gene.symbol,
        organism="human",
    )
    validated = curator.validate()
    assert validated is False
    assert curator.non_validated == {
        "cell_type": ["astrocytic glia"],
        "donor": ["D0001", "D0002", "DOOO3"],
        "var_index": ["TCF-1"],
    }

    # standardize var_index
    curator.standardize("var_index")
    assert "TCF7" in adata.var.index
    assert curator.non_validated == {
        "cell_type": ["astrocytic glia"],
        "donor": ["D0001", "D0002", "DOOO3"],
    }
    curator.standardize("all")
    assert curator.non_validated == {"donor": ["D0001", "D0002", "DOOO3"]}

    # lookup
    lookup = curator.lookup()
    assert lookup.cell_type.oligodendrocyte.name == "oligodendrocyte"

    # add new
    curator.add_new_from(to_add)
    assert curator.non_validated == {}
    # just for coverage, doesn't do anything
    curator.add_new_from_var_index()
    validated = curator.validate()
    assert validated

    artifact = curator.save_artifact(description="test AnnData")

    assert set(artifact.features.get_values()["cell_type"]) == {
        "cerebral cortex pyramidal neuron",
        "astrocyte",
        "oligodendrocyte",
    }
    assert set(artifact.features.get_values()["cell_type_2"]) == {
        "oligodendrocyte",
        "astrocyte",
    }

    # clean up
    artifact.delete(permanent=True)
    ln.ULabel.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.CellType.filter().delete()
    ln.FeatureSet.filter().delete()
    bt.Gene.filter().delete()


def test_str_var_index(adata):
    with pytest.raises(TypeError, match="var_index parameter has to be a bionty field"):
        _ = ln.Curator.from_anndata(
            adata,
            var_index="symbol",
            organism="human",
        )


def test_not_passing_categoricals(adata):
    curator = ln.Curator.from_anndata(
        adata,
        var_index=bt.Gene.symbol,
        organism="human",
    )
    validated = curator.validate()
    assert validated is False


def test_anndata_curator_wrong_type(df, categoricals):
    with pytest.raises(TypeError, match="data has to be an AnnData object"):
        ln.Curator.from_anndata(
            df,
            categoricals=categoricals,
            var_index=bt.Gene.symbol,
            organism="human",
        )


def test_categorical_key_not_present(df):
    with pytest.raises(
        ValidationError,
        match="the following 1 key passed to categoricals is not allowed:",
    ):
        ln.Curator.from_df(
            df,
            categoricals={"not present": None},
            organism="human",
        )


def test_source_key_not_present(adata, categoricals):
    with pytest.raises(
        ValidationError,
        match="the following 1 key passed to sources is not allowed:",
    ):
        ln.Curator.from_anndata(
            adata,
            categoricals=categoricals,
            var_index=bt.Gene.symbol,
            sources={"not_present": None},
            organism="human",
        )


def test_unvalidated_adata_object(adata, categoricals):
    curator = ln.Curator.from_anndata(
        adata,
        categoricals=categoricals,
        var_index=bt.Gene.symbol,
        organism="human",
    )
    with pytest.raises(
        ValidationError, match="Dataset does not validate. Please curate."
    ):
        curator.save_artifact()


def test_mudata_curator(mdata):
    mdata = mdata.copy()
    categoricals = {
        "rna:cell_type": bt.CellType.name,
        "rna:assay_ontology_id": bt.ExperimentalFactor.ontology_id,
        "rna:donor": ln.ULabel.name,
        "rna_2:cell_type": bt.CellType.name,
        "rna_2:assay_ontology_id": bt.ExperimentalFactor.ontology_id,
        "rna_2:donor": ln.ULabel.name,
        "donor": ln.ULabel.name,
    }

    curator = ln.Curator.from_mudata(
        mdata,
        categoricals=categoricals,
        var_index={"rna": bt.Gene.symbol, "rna_2": bt.Gene.symbol},
        organism="human",
    )
    with pytest.raises(ValidationError):
        _ = curator.non_validated
    assert curator._modalities == curator._modalities

    # validate
    validated = curator.validate()
    assert curator.non_validated == {
        "obs": {"donor": ["D0001", "D0002", "DOOO3"]},
        "rna_2": {
            "cell_type": ["astrocytic glia"],
            "donor": ["D0001", "D0002", "DOOO3"],
            "var_index": ["TCF-1"],
        },
        "rna": {
            "cell_type": ["astrocytic glia"],
            "donor": ["D0001", "D0002", "DOOO3"],
            "var_index": ["TCF-1"],
        },
    }

    # lookup
    _ = curator.lookup()

    # standardize
    curator.standardize("all", modality="rna")
    curator.standardize("all", modality="rna_2")
    assert curator._mod_adata_curators["rna_2"].non_validated == {
        "donor": ["D0001", "D0002", "DOOO3"]
    }

    # add new
    curator.add_new_from_var_index("rna")  # doesn't do anything
    curator.add_new_from("donor")

    validated = curator.validate()
    assert validated
    artifact = curator.save_artifact(description="test MuData")

    # clean up
    artifact.delete(permanent=True)
    ln.ULabel.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.CellType.filter().delete()
    ln.FeatureSet.filter().delete()
    bt.Gene.filter().delete()
