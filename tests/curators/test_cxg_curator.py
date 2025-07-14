import bionty as bt
import lamindb as ln
import pytest


def test_cxg_curator():
    ln.examples.cellxgene.save_cxg_defaults()

    schema = ln.examples.cellxgene.get_cxg_schema(
        schema_version="5.2.0", field_types=["name", "ontology_id"]
    )
    adata = ln.core.datasets.small_dataset3_cellxgene(with_defaults=True)

    curator = ln.curators.AnnDataCurator(adata, schema)
    # Ensure that default values for Features are set
    curator.slots["obs"].standardize()

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
    features = [f for sm in schema.slots.values() for f in sm.features.all()]
    ln.models.SchemaFeature.filter(feature__in=features).delete()
    schema.delete()
    [f.delete() for f in features]

    for ulabel in ln.ULabel.filter(is_type=False):
        ulabel.delete()
    for entity in [
        bt.Disease,
        bt.Ethnicity,
        bt.DevelopmentalStage,
        bt.Phenotype,
        bt.CellType,
        ln.ULabel,
    ]:
        for record in entity.filter().all():
            record.delete()


def test_invalid_field_type():
    with pytest.raises(ValueError) as e:
        ln.examples.cellxgene.get_cxg_schema(
            schema_version="5.3.0", field_types=["ensembl_gene_ids"]
        )
    assert "Invalid field_types" in str(e.value)


def test_invalid_schema_Version():
    with pytest.raises(ValueError) as e:
        ln.examples.cellxgene.get_cxg_schema(schema_version="200.0.0")
    assert "Invalid schema_version" in str(e.value)
