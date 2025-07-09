import bionty as bt
import lamindb as ln
import pytest


def test_cxg_curator():
    ln.examples.cellxgene.save_cxg_defaults()

    schema = ln.examples.cellxgene.get_cxg_schema(
        schema_version="5.2.0", key_types=["name", "ontology_id"]
    )
    adata = ln.core.datasets.small_dataset3_cellxgene()

    curator = ln.curators.AnnDataCurator(adata, schema)

    ln.examples.cellxgene.add_cxg_defaults_to_obs(adata)

    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()
        assert (
            "ValidationError: 1 term not validated in feature 'index' in slot 'var': 'invalid_ensembl_id'"
            in str(e)
        )

    adata = adata[
        :, ~adata.var.index.isin(curator.slots["var"].cat.non_validated["index"])
    ].copy()

    curator = ln.curators.AnnDataCurator(adata, schema)

    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()
        assert (
            "ValidationError: 1 term not validated in feature 'tissue' in slot 'obs': 'lungg'"
            in str(e)
        )

    adata.obs["tissue"] = adata.obs["tissue"].cat.rename_categories({"lungg": "lung"})

    curator.validate()

    artifact = curator.save_artifact(key="examples/dataset-curated-against-cxg.h5ad")

    # Clean up
    artifact.delete(permanent=True)
    bt.Disease.get(uid=bt.Phenotype.get(ontology_id="PATO:0000461").uid).delete()
    for model, name in zip(
        [
            bt.Ethnicity,
            bt.Ethnicity,
            bt.DevelopmentalStage,
            bt.Phenotype,
            bt.CellType,
        ],
        ["na", "unknown", "unknown", "unknown", "unknown"],
    ):
        model.get(
            ontology_id=name, name=name, description="From CellxGene schema."
        ).delete()
    for name in ["tissue", "organoid", "cell culture"]:
        ln.ULabel.get(name=name, description="From CellxGene schema.").delete()
    ln.ULabel.get(name="TissueType", is_type=True).delete()
    for name in ["cell", "nucleus", "na"]:
        ln.ULabel.get(name=name, description="From CellxGene schema.").delete()
    ln.ULabel.get(
        name="SuspensionType",
        is_type=True,
    ).delete()
