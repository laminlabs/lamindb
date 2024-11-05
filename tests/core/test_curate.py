from unittest.mock import Mock

import anndata as ad
import bionty as bt
import lamindb as ln
import mudata as md
import pandas as pd
import pytest
from lamindb._curate import CurateLookup
from lamindb.core.exceptions import ValidationError


@pytest.fixture(scope="module")
def df():
    return pd.DataFrame(
        {
            "cell_type": [
                # there is an error in the below annotation on purpose
                "cerebral pyramidal neuron",
                "astrocyte",
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
                "astrocyte",
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
            "TCF7": [1, 2, 3],
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
    curate = ln.Curator.from_df(df, categoricals=categoricals)
    validated = curate.validate()
    assert validated is False

    cell_types = curate.lookup(public=True)["cell_type"]
    df["cell_type"] = df["cell_type"].replace(
        {"cerebral pyramidal neuron": cell_types.cerebral_cortex_pyramidal_neuron.name}
    )
    curate.add_new_from("donor")
    validated = curate.validate()
    assert validated is True

    artifact = curate.save_artifact(description="test-curate-df")

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


def test_custom_using_invalid_field_lookup(curate_lookup):
    with pytest.raises(AttributeError) as excinfo:
        _ = curate_lookup["invalid_field"]
    assert "'CurateLookup' object has no attribute 'invalid_field'" in str(
        excinfo.value
    )


def test_additional_args_with_all_key(df, categoricals):
    curate = ln.Curator.from_df(df, categoricals=categoricals)
    with pytest.raises(ValueError) as error:
        curate.add_new_from("all", extra_arg="not_allowed")
    assert "Cannot pass additional arguments to 'all' key!" in str(error.value)


def test_save_columns_not_defined_in_fields(df, categoricals):
    curate = ln.Curator.from_df(df, categoricals=categoricals)
    with pytest.raises(ValidationError) as error:
        curate._update_registry("nonexistent")
    assert "Feature nonexistent is not part of the fields!" in str(error.value)


def test_unvalidated_data_object(df, categoricals):
    curate = ln.Curator.from_df(df, categoricals=categoricals)
    with pytest.raises(ValidationError) as error:
        curate.save_artifact()
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

    curate = ln.Curator.from_df(pd.DataFrame())
    curate.clean_up_failed_runs()

    assert len(ln.Run.filter(transform=mock_transform).all()) == 1

    # Revert to old run context to not infer with tests that need the run context
    context._transform = previous_transform
    context._run = previous_run


@pytest.mark.parametrize("to_add", ["donor", "all"])
def test_anndata_curator(adata, categoricals, to_add):
    curate = ln.Curator.from_anndata(
        adata,
        categoricals=categoricals,
        var_index=bt.Gene.symbol,
        organism="human",
    )
    curate.add_new_from(to_add)
    validated = curate.validate()
    assert validated

    artifact = curate.save_artifact(description="test AnnData")

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


def test_no_categoricals(adata):
    curate = ln.Curate.from_anndata(
        adata,
        var_index=bt.Gene.symbol,
        organism="human",
    )
    validated = curate.validate()
    assert validated


def test_anndata_curator_wrong_type(df, categoricals):
    with pytest.raises(ValueError) as error:
        ln.Curator.from_anndata(
            df,
            categoricals=categoricals,
            var_index=bt.Gene.symbol,
            organism="human",
        )
    assert "data has to be an AnnData object" in str(error.value)


def test_categorical_key_not_present(df):
    with pytest.raises(ValidationError) as error:
        ln.Curator.from_df(
            df,
            categoricals={"not present": None},
            organism="human",
        )

    assert "the following keys passed to categoricals are not allowed:" in str(
        error.value
    )


def test_source_key_not_present(adata, categoricals):
    with pytest.raises(ValidationError) as error:
        ln.Curator.from_anndata(
            adata,
            categoricals=categoricals,
            var_index=bt.Gene.symbol,
            sources={"not_present": None},
            organism="human",
        )
    assert "the following keys passed to sources are not allowed:" in str(error.value)


def test_unvalidated_adata_object(adata, categoricals):
    curate = ln.Curator.from_anndata(
        adata,
        categoricals=categoricals,
        var_index=bt.Gene.symbol,
        organism="human",
    )
    with pytest.raises(ValidationError) as error:
        curate.save_artifact()
    assert "Dataset does not validate. Please curate." in str(error.value)


def test_mudata_curator(mdata):
    categoricals = {
        "rna:cell_type": bt.CellType.name,
        "rna:assay_ontology_id": bt.ExperimentalFactor.ontology_id,
        "rna:donor": ln.ULabel.name,
        "rna_2:cell_type": bt.CellType.name,
        "rna_2:assay_ontology_id": bt.ExperimentalFactor.ontology_id,
        "rna_2:donor": ln.ULabel.name,
    }

    curate = ln.Curator.from_mudata(
        mdata,
        categoricals=categoricals,
        var_index={"rna": bt.Gene.symbol, "rna_2": bt.Gene.symbol},
        organism="human",
    )
    curate.add_new_from("donor", modality="rna")
    validated = curate.validate()
    assert validated
    artifact = curate.save_artifact(description="test MuData")

    # clean up
    artifact.delete(permanent=True)
    ln.ULabel.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.CellType.filter().delete()
