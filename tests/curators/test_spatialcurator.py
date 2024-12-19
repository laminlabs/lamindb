import re

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from spatialdata.datasets import blobs


@pytest.fixture
def blobs_data():
    sdata = blobs()

    sdata.tables["table"].var.index = [
        "ENSG00000139618",  # BRCA2
        "ENSG00000157764",  # BRAF
        "ENSG00000999999",  # does not exist - to test add_new_from_var_index
    ]
    sdata.tables["table"].obs["region"] = pd.Categorical(
        ["region 1"] * 13 + ["region 2"] * 13
    )
    sdata.attrs["sample"] = {
        "assay": "Visium Spatial Gene Expression",
        "disease": "Alzheimer's dementia",
        "developmental_stage": "very early",  # does not exist - to test add_new_from
    }

    return sdata


def test_spatialdata_curator(blobs_data):
    from lamin_spatial import SpatialDataCurator
    from lamindb.core.exceptions import ValidationError

    with pytest.raises(
        ValidationError, match="key passed to categoricals is not present"
    ):
        SpatialDataCurator(
            blobs_data,
            var_index={"table": bt.Gene.ensembl_gene_id},
            categoricals={
                "sample": {
                    "does not exist": bt.ExperimentalFactor.name,
                },
            },
            organism="human",
        )

    with pytest.raises(ValidationError, match="key passed to sources is not present"):
        SpatialDataCurator(
            blobs_data,
            var_index={"table": bt.Gene.ensembl_gene_id},
            categoricals={
                "table": {"region": ln.ULabel.name},
            },
            sources={"sample": {"whatever": bt.CellLine.name}},
            organism="human",
        )

    curator = SpatialDataCurator(
        blobs_data,
        var_index={"table": bt.Gene.ensembl_gene_id},
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

    with pytest.raises(ValidationError, match=re.escape("Run .validate() first.")):
        curator.add_new_from(key="region", accessor="table")

    with pytest.raises(
        ValidationError, match=re.escape("Dataset does not validate. Please curate.")
    ):
        curator.save_artifact(description="test spatialdata curation")

    assert not curator.validate()

    assert curator.non_validated == {
        "sample": {
            "disease": ["Alzheimer's dementia"],
            "developmental_stage": ["very early"],
        },
        "table": {"region": ["region 1", "region 2"], "var_index": ["ENSG00000999999"]},
    }

    curator.add_new_from_var_index("table")

    assert curator.non_validated == {
        "sample": {
            "disease": ["Alzheimer's dementia"],
            "developmental_stage": ["very early"],
        },
        "table": {"region": ["region 1", "region 2"]},
    }

    curator.add_new_from(key="developmental_stage", accessor="sample")
    curator.add_new_from(key="region", accessor="table")

    assert curator.non_validated == {"sample": {"disease": ["Alzheimer's dementia"]}}

    # test invalid key in standardize
    with pytest.raises(ValueError, match="key 'invalid_key' not present in 'table'!"):
        curator.standardize(key="invalid_key", accessor="table")

    # standardize
    assert curator.non_validated == {"sample": {"disease": ["Alzheimer's dementia"]}}
    curator.standardize(key="disease", accessor="sample")
    assert curator._sample_metadata["disease"].values[0] == "Alzheimer disease"
    assert curator.non_validated == {}

    # validation should finally pass
    assert curator.validate() is True

    # lookup
    lookup = curator.lookup()
    assert lookup.disease[0].name == "Alzheimer disease"

    # save & associated features
    artifact = curator.save_artifact(description="blob spatialdata")
    assert artifact.features.get_values()["assay"] == "Visium Spatial Gene Expression"
    assert set(artifact.features.get_values()["region"]) == {"region 1", "region 2"}

    # clean up
    artifact.delete(permanent=True)
    ln.ULabel.filter().delete()
    bt.ExperimentalFactor.filter().delete()
    bt.Disease.filter().delete()
    bt.DevelopmentalStage.filter().delete()
    ln.FeatureSet.filter().delete()
    bt.Gene.filter().delete()
    ln.Feature.filter().delete()
