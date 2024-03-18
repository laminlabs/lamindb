import anndata as ad
import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.validation import AnnDataValidator, Validator


@pytest.fixture(scope="module")
def df():
    return pd.DataFrame(
        {
            "cell_type": ["cerebral pyramidal neuron", "astrocyte", "oligodendrocyte"],
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
def fields():
    return {
        "cell_type": bt.CellType.name,
        "assay_ontology_id": bt.ExperimentalFactor.ontology_id,
        "donor": ln.ULabel.name,
    }


def test_validator(df, fields):
    validator = Validator(df, fields=fields)
    validated = validator.validate()
    assert validated is False

    cell_types = validator.lookup("public")["cell_type"]
    df["cell_type"] = df["cell_type"].replace(
        {"cerebral pyramidal neuron": cell_types.cerebral_cortex_pyramidal_neuron.name}
    )
    validator.register_labels("all")
    validator.register_labels("donor", validated_only=False)
    validated = validator.validate()
    assert validated is True


def test_anndata_validator(adata, fields):
    validator = AnnDataValidator(
        adata,
        obs_fields=fields,
        var_field=bt.Gene.symbol,  # specify the field for the var
    )
    validated = validator.validate(organism="human")
    assert validated is False

    validator.register_labels("variables")
    validated = validator.validate()
    assert validated is True

    artifact = validator.register_artifact(description="test AnnData")
    collection = validator.register_collection(
        artifact,
        name="Experiment X in brain",
        description="10.1126/science.xxxxx",
        reference="E-MTAB-xxxxx",
        reference_type="ArrayExpress",
    )
    assert collection.artifact == artifact
