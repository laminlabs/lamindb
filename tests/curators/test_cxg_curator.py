import bionty as bt
import lamindb as ln
import pytest


@pytest.fixture(scope="module", params=["5.2.0", "6.0.0"])
def cxg_schema(request):
    ln.examples.cellxgene.save_cxg_defaults()
    schema = ln.examples.cellxgene.get_cxg_schema(
        schema_version=request.param, field_types=["name", "ontology_id"]
    )

    yield schema, request.param

    ln.models.SchemaComponent.filter().delete()
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    ln.ULabel.filter(type__isnull=False).delete()
    for entity in [
        bt.Disease,
        bt.Ethnicity,
        bt.DevelopmentalStage,
        bt.Phenotype,
        bt.CellType,
        ln.ULabel,
    ]:
        entity.filter().all().delete()


def test_cxg_curator(cxg_schema):
    schema, version = cxg_schema
    is_v6 = version == "6.0.0"

    # test invalid var index and typo in obs column
    adata = ln.core.datasets.small_dataset3_cellxgene(
        with_obs_defaults=True, with_obs_typo=True, with_uns_organism=is_v6
    )
    curator = ln.curators.AnnDataCurator(adata, schema)
    # Ensure that default values for Features are set
    curator.slots["obs"].standardize()
    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()
    assert (
        "1 term not validated in feature 'index' in slot 'var': 'invalid_ensembl_id'"
        in str(e)
    )
    # subset adata to remove invalid var index
    adata = adata[
        :, ~adata.var.index.isin(curator.slots["var"].cat.non_validated["index"])
    ].copy()
    curator = ln.curators.AnnDataCurator(adata, schema)
    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()
    assert (
        "1 term not validated in feature 'tissue_ontology_term_id' in slot 'obs': 'UBERON:0002048XXX'"
        in str(e)
    )
    # fix typo in obs column
    adata.obs["tissue_ontology_term_id"] = adata.obs["tissue_ontology_term_id"].replace(
        {"UBERON:0002048XXX": "UBERON:0002048"}
    )
    artifact = curator.save_artifact(
        key=f"examples/dataset-curated-against-cxg-{version}.h5ad"
    )

    # test missing obs columns
    adata = ln.core.datasets.small_dataset3_cellxgene(
        with_obs_defaults=False, with_uns_organism=is_v6
    )
    adata = adata[:, ~adata.var.index.isin({"invalid_ensembl_id"})].copy()
    curator = ln.curators.AnnDataCurator(adata, schema)
    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()

    expected_missing_obs = [
        "assay",
        "cell_type",
        "development_stage",
        "disease",
        "self_reported_ethnicity",
        "assay_ontology_term_id",
        "cell_type_ontology_term_id",
        "development_stage_ontology_term_id",
        "self_reported_ethnicity_ontology_term_id",
        "tissue_ontology_term_id",
        "donor_id",
    ]
    if not is_v6:
        expected_missing_obs.append("organism_ontology_term_id")

    assert all(col in str(e.value) for col in expected_missing_obs)

    # Clean up
    artifact.delete(permanent=True)


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
