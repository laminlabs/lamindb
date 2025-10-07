import re
from unittest.mock import Mock

import anndata as ad
import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.curators.core import CatLookup, ValidationError


@pytest.fixture
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
            "donor": ["D0001", "D0002", "D0003"],
        }
    )


@pytest.fixture
def adata():
    # this should be using small_dataset1 instead of the custom code here
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
            "donor": ["D0001", "D0002", "D0003"],
            "sample_note": ["was ok", "looks naah", "pretty! 🤩"],
            "temperature": [23.1, 23.2, 33.3],
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


@pytest.fixture
def mdata(adata):
    import mudata as md

    # can't be the same adata object due to in-place modifications
    mdata = md.MuData({"rna": adata, "rna_2": adata.copy()})
    mdata.obs["donor"] = ["D0001", "D0002", "D0003"]

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
    return CatLookup(categoricals=categoricals)


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


def test_df_cat_manager(df):
    try:
        cat_manager = ln.curators.core.DataFrameCatManager(
            df,
            categoricals=[
                ln.Feature(name="cell_type", dtype=bt.CellType).save(),
                ln.Feature(name="cell_type_2", dtype=bt.CellType).save(),
                ln.Feature(
                    name="assay_ontology_id", dtype=bt.ExperimentalFactor.ontology_id
                ).save(),
                ln.Feature(name="donor", dtype=ln.ULabel).save(),
            ],
        )
        with pytest.raises(ValidationError):
            _ = cat_manager.non_validated
        validated = cat_manager.validate()
        assert cat_manager.non_validated == {
            "cell_type": ["cerebral pyramidal neuron", "astrocytic glia"],
            "donor": ["D0001", "D0002", "D0003"],
        }
        assert validated is False

        # standardize
        with pytest.raises(KeyError):
            cat_manager.standardize("nonexistent-key")
        cat_manager.standardize("all")
        assert cat_manager.non_validated == {
            "cell_type": ["cerebral pyramidal neuron"],
            "donor": ["D0001", "D0002", "D0003"],
        }
        assert "astrocyte" in df["cell_type"].values

        # add new
        cat_manager.add_new_from("donor")
        assert cat_manager.non_validated == {"cell_type": ["cerebral pyramidal neuron"]}

        # lookup
        cell_types = cat_manager.lookup(public=True)["cell_type"]
        df["cell_type"] = df["cell_type"].replace(
            {
                "cerebral pyramidal neuron": cell_types.cerebral_cortex_pyramidal_neuron.name
            }
        )
        validated = cat_manager.validate()
        assert validated is True
        assert cat_manager.non_validated == {}

        # no need to standardize
        cat_manager.standardize("cell_type")
    finally:
        ln.ULabel.filter().delete(permanent=True)
        bt.ExperimentalFactor.filter().delete(permanent=True)
        bt.CellType.filter().delete(permanent=True)
        ln.Schema.filter().delete(permanent=True)


def test_custom_using_invalid_field_lookup(curate_lookup):
    with pytest.raises(
        AttributeError, match='"CatLookup" object has no attribute "invalid_field"'
    ):
        _ = curate_lookup["invalid_field"]


def test_additional_args_with_all_key(df, categoricals):
    curator = ln.Curator.from_dataframe(df, categoricals=categoricals)
    with pytest.raises(
        ValueError, match="Cannot pass additional arguments to 'all' key!"
    ):
        curator.add_new_from("all", extra_arg="not_allowed")


# below is all deprecated


def test_unvalidated_data_object(df, categoricals):
    curator = ln.Curator.from_dataframe(df, categoricals=categoricals)
    with pytest.raises(
        ValidationError, match="Dataset does not validate. Please curate."
    ):
        curator.save_artifact()


def test_df_curator(df, categoricals):
    try:
        curator = ln.Curator.from_dataframe(df, categoricals=categoricals)
        with pytest.raises(ValidationError):
            _ = curator.non_validated
        validated = curator.validate()
        assert curator.non_validated == {
            "cell_type": ["cerebral pyramidal neuron", "astrocytic glia"],
            "donor": ["D0001", "D0002", "D0003"],
        }
        assert validated is False

        # standardize
        with pytest.raises(KeyError):
            curator.standardize("nonexistent-key")
        curator.standardize("all")
        assert curator.non_validated == {
            "cell_type": ["cerebral pyramidal neuron"],
            "donor": ["D0001", "D0002", "D0003"],
        }
        assert "astrocyte" in df["cell_type"].values

        # add new
        curator.add_new_from("donor")
        assert curator.non_validated == {"cell_type": ["cerebral pyramidal neuron"]}

        # lookup
        cell_types = curator.lookup(public=True)["cell_type"]
        df["cell_type"] = df["cell_type"].replace(
            {
                "cerebral pyramidal neuron": cell_types.cerebral_cortex_pyramidal_neuron.name
            }
        )
        validated = curator.validate()
        assert validated is True
        assert curator.non_validated == {}

        # no need to standardize
        curator.standardize("cell_type")

        artifact = curator.save_artifact(description="test-curate-df")

        artifact.describe()

        assert (
            artifact.cell_types.through.filter(artifact_id=artifact.id)
            .to_dataframe()["label_ref_is_name"]
            .values.sum()
            == 5
        )

        assert (
            artifact.experimental_factors.through.filter(artifact_id=artifact.id)
            .to_dataframe()["label_ref_is_name"]
            .values.sum()
            == 0
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
    finally:
        # clean up
        artifact.delete(permanent=True)
        ln.ULabel.filter().delete(permanent=True)
        bt.ExperimentalFactor.filter().delete(permanent=True)
        bt.CellType.filter().delete(permanent=True)
        ln.Schema.filter().delete(permanent=True)


def test_pass_artifact(df):
    try:
        artifact = ln.Artifact.from_dataframe(
            df, key="test_cat_curators/df.parquet"
        ).save()
        curator = ln.Curator.from_dataframe(
            artifact, categoricals={"donor": ln.ULabel.name}
        )
        curator.validate()
        with pytest.raises(
            RuntimeError,
            match="can't mutate the dataset when an artifact is passed!",
        ):
            curator.standardize("all")
        curator.add_new_from("donor")
        artifact_2 = curator.save_artifact()
        assert artifact == artifact_2
    finally:
        # clean up
        artifact.delete(permanent=True)
        ln.ULabel.filter().delete(permanent=True)
        ln.Schema.filter().delete(permanent=True)


@pytest.mark.parametrize("to_add", ["donor", "all"])
def test_anndata_curator(adata, categoricals, to_add):
    try:
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
            "donor": ["D0001", "D0002", "D0003"],
            "var_index": ["TCF-1"],
        }

        # standardize var_index
        curator.standardize("var_index")

        assert "TCF7" in adata.var.index
        assert curator.non_validated == {
            "cell_type": ["astrocytic glia"],
            "donor": ["D0001", "D0002", "D0003"],
        }
        # deprecated
        curator.standardize("all")
        assert curator.non_validated == {"donor": ["D0001", "D0002", "D0003"]}

        # lookup
        lookup = curator.lookup()
        assert lookup.cell_type.oligodendrocyte.name == "oligodendrocyte"

        # add new
        curator.add_new_from(to_add)
        assert curator.non_validated == {}
        # deprecated
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
    finally:
        # clean up
        artifact.delete(permanent=True)
        ln.ULabel.filter().delete(permanent=True)
        bt.ExperimentalFactor.filter().delete(permanent=True)
        bt.CellType.filter().delete(permanent=True)
        ln.Schema.filter().delete(permanent=True)
        bt.Gene.filter().delete(permanent=True)


def test_str_var_index(adata):
    with pytest.raises(
        TypeError,
        match="var_index parameter has to be a field, e.g. Gene.ensembl_gene_id",
    ):
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
    categoricals = {
        "rna:cell_type": bt.CellType.name,
        "rna:assay_ontology_id": bt.ExperimentalFactor.ontology_id,
        "rna:donor": ln.ULabel.name,
        "rna_2:cell_type": bt.CellType.name,
        "rna_2:assay_ontology_id": bt.ExperimentalFactor.ontology_id,
        "rna_2:donor": ln.ULabel.name,
        "donor": ln.ULabel.name,
    }

    try:
        artifact = None
        curator = ln.Curator.from_mudata(
            mdata,
            categoricals=categoricals,
            var_index={"rna": bt.Gene.symbol, "rna_2": bt.Gene.symbol},
            organism="human",
        )
        with pytest.raises(ValidationError):
            _ = curator.non_validated
        assert curator._modalities == {"obs", "rna", "rna_2"}

        # validate
        validated = curator.validate()
        assert curator.non_validated == {
            "obs": {"donor": ["D0001", "D0002", "D0003"]},
            "rna_2": {
                "cell_type": ["astrocytic glia"],
                "donor": ["D0001", "D0002", "D0003"],
                "var_index": ["TCF-1"],
            },
            "rna": {
                "cell_type": ["astrocytic glia"],
                "donor": ["D0001", "D0002", "D0003"],
                "var_index": ["TCF-1"],
            },
        }

        # lookup
        curator.lookup()

        # standardize
        curator.standardize("all", modality="rna")
        curator.standardize("all", modality="rna_2")
        assert curator._mod_adata_curators["rna_2"].non_validated == {
            "donor": ["D0001", "D0002", "D0003"]
        }

        # add new
        curator.add_new_from_var_index("rna")  # doesn't do anything
        curator.add_new_from("donor")

        validated = curator.validate()
        assert validated
        artifact = curator.save_artifact(description="test MuData")
    finally:
        # clean up
        if artifact:
            artifact.delete(permanent=True)
        ln.ULabel.filter().delete(permanent=True)
        bt.ExperimentalFactor.filter().delete(permanent=True)
        bt.CellType.filter().delete(permanent=True)
        ln.Schema.filter().delete(permanent=True)
        bt.Gene.filter().delete(permanent=True)


def test_spatialdata_curator():
    from spatialdata.datasets import blobs

    blobs_data = blobs()

    blobs_data.tables["table"].var.index = [
        "TSPAN6",
        "MYODULIN",  # synonym
        "DOESNOTEXIST",  # does not exist - to test add_new_from_var_index
    ]
    blobs_data.tables["table"].obs["region"] = pd.Categorical(
        ["region 1"] * 13 + ["region 2"] * 13
    )
    blobs_data.attrs["sample"] = {
        "assay": "Visium Spatial Gene Expression",
        "disease": "Alzheimer's dementia",
        "developmental_stage": "very early",  # does not exist - to test add_new_from
    }

    from lamindb.errors import ValidationError

    with pytest.raises(
        ValidationError, match="key passed to categoricals is not present"
    ):
        ln.Curator.from_spatialdata(
            blobs_data,
            var_index={"table": bt.Gene.symbol},
            categoricals={
                "sample": {
                    "does not exist": bt.ExperimentalFactor.name,
                },
            },
            organism="human",
        )

    with pytest.raises(ValidationError, match="key passed to sources is not present"):
        ln.Curator.from_spatialdata(
            blobs_data,
            var_index={"table": bt.Gene.symbol},
            categoricals={
                "table": {"region": ln.ULabel.name},
            },
            sources={"sample": {"whatever": bt.CellLine.name}},
            organism="human",
        )

    try:
        curator = ln.Curator.from_spatialdata(
            blobs_data,
            var_index={"table": bt.Gene.symbol},
            categoricals={
                "sample": {
                    "assay": bt.ExperimentalFactor.name,
                    "disease": bt.Disease.name,
                    "developmental_stage": bt.DevelopmentalStage.name,
                },
                "table": {"region": ln.ULabel.name},
            },
            organism="human",
        )

        with pytest.raises(
            ValidationError, match=re.escape("Please run validate() first!")
        ):
            curator.add_new_from(key="region", accessor="table")

        with pytest.raises(
            ValidationError,
            match=re.escape("Dataset does not validate. Please curate."),
        ):
            curator.save_artifact(description="test spatialdata curation")

        assert not curator.validate()

        assert curator.non_validated == {
            "sample": {
                "disease": ["Alzheimer's dementia"],
                "developmental_stage": ["very early"],
            },
            "table": {
                "region": ["region 1", "region 2"],
                "var_index": ["MYODULIN", "DOESNOTEXIST"],
            },
        }

        curator.add_new_from(key="developmental_stage", accessor="sample")
        curator.add_new_from(key="region", accessor="table")

        assert curator.non_validated == {
            "sample": {"disease": ["Alzheimer's dementia"]},
            "table": {"var_index": ["MYODULIN", "DOESNOTEXIST"]},
        }

        # test invalid key in standardize
        with pytest.raises(
            ValueError, match="key 'invalid_key' not present in 'table'!"
        ):
            curator.standardize(key="invalid_key", accessor="table")

        # standardize
        curator.standardize(key="disease", accessor="sample")
        assert curator._sample_metadata["disease"].values[0] == "Alzheimer disease"
        curator.standardize(key="var_index", accessor="table")
        assert curator.non_validated == {"table": {"var_index": ["DOESNOTEXIST"]}}
        curator.add_new_from(key="var_index", accessor="table")
        assert curator.non_validated == {}

        # validation should finally pass
        assert curator.validate() is True

        # lookup
        lookup = curator.lookup()
        assert lookup.disease[0].name == "Alzheimer disease"

        # save & associated features
        artifact = curator.save_artifact(description="blob spatialdata")
        # the two below tests broke in https://github.com/laminlabs/lamindb/pull/2650
        # but only for the legacy curator
        # assert (
        #     artifact.features.get_values()["assay"] == "Visium Spatial Gene Expression"
        # )
        # assert set(artifact.features.get_values()["region"]) == {"region 1", "region 2"}

    finally:
        artifact.delete(permanent=True)
        ln.ULabel.filter().delete(permanent=True)
        bt.ExperimentalFactor.filter().delete(permanent=True)
        bt.Disease.filter().delete(permanent=True)
        bt.DevelopmentalStage.filter().delete(permanent=True)
        ln.Schema.filter().delete(permanent=True)
        bt.Gene.filter().delete(permanent=True)
        ln.Feature.filter().delete(permanent=True)
