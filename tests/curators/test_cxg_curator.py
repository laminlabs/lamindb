import bionty as bt
import lamindb as ln
import pytest
from lamindb.examples.cellxgene._cellxgene import CELLxGENESchemaVersions


@pytest.fixture
def cxg_schema_factory():
    def create_schema(version: CELLxGENESchemaVersions, **kwargs):
        ln.examples.cellxgene.save_cellxgene_defaults()
        schema = ln.examples.cellxgene.create_cellxgene_schema(
            schema_version=version, **kwargs
        )
        return schema

    yield create_schema

    # Cleanup after all tests
    ln.models.SchemaComponent.filter().delete(permanent=True)
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    ln.ULabel.filter(type__isnull=False).delete(permanent=True)
    for entity in [
        bt.Disease,
        bt.Ethnicity,
        bt.DevelopmentalStage,
        bt.Phenotype,
        bt.CellType,
        ln.ULabel,
    ]:
        entity.filter().delete(permanent=True)


def test_cxg_curator_5(cxg_schema_factory):
    cxg_schema = cxg_schema_factory("5.2.0", field_types=["name", "ontology_id"])

    # test invalid var index and typo in obs column
    adata = ln.examples.datasets.small_dataset3_cellxgene(
        with_obs_defaults=True, with_obs_typo=True, with_var_typo=True
    )
    curator = ln.curators.AnnDataCurator(adata, cxg_schema)
    # Ensure that default values for Features are set
    curator.slots["obs"].standardize()
    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()
        assert (
            "ValidationError: 1 term not validated in feature 'index' in slot 'var': 'invalid_ensembl_id'"
            in str(e)
        )

    # subset adata to remove invalid var index
    adata = adata[
        :, ~adata.var.index.isin(curator.slots["var"].cat.non_validated["index"])
    ].copy()
    curator = ln.curators.AnnDataCurator(adata, cxg_schema)
    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()
        assert (
            "ValidationError: 1 term not validated in feature 'tissue_ontology_term_id' in slot 'obs': 'UBERON:0002048XXX'"
            in str(e)
        )
    # fix typo in obs column
    adata.obs["tissue_ontology_term_id"] = adata.obs["tissue_ontology_term_id"].replace(
        {"UBERON:0002048XXX": "UBERON:0002048"}
    )
    artifact = curator.save_artifact(key="examples/dataset-curated-against-cxg-5.h5ad")

    # test missing obs columns
    adata = ln.examples.datasets.small_dataset3_cellxgene(with_obs_defaults=False)
    adata = adata[:, ~adata.var.index.isin({"invalid_ensembl_id"})].copy()
    curator = ln.curators.AnnDataCurator(adata, cxg_schema)
    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()
    expected_missing = [
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
        "organism_ontology_term_id",
        "donor_id",
    ]
    assert all(col in str(e.value) for col in expected_missing)

    # Clean up
    artifact.delete(permanent=True)


def test_cxg_curator_6_spatial(cxg_schema_factory):
    """Tests organism (in `uns` as of 6.x) and spatial slot validation of CELLxGENE 6.x."""
    cxg_schema = cxg_schema_factory(
        "6.0.0", spatial_library_id="library_123", field_types="ontology_id"
    )

    adata = ln.examples.datasets.small_dataset3_cellxgene(
        with_obs_defaults=True, with_uns_organism=True, with_uns_spatial=True
    )
    # delete a necessary component from uns["spatial"]
    del adata.uns["spatial"]["is_single"]

    curator = ln.curators.AnnDataCurator(adata, cxg_schema)

    with pytest.raises(ln.errors.ValidationError) as e:
        curator.validate()
    assert "column 'is_single' not in dataframe." in str(e.value)

    adata.uns["spatial"]["is_single"] = True
    curator = ln.curators.AnnDataCurator(adata, cxg_schema)
    curator.validate()


def test_invalid_field_type():
    with pytest.raises(ValueError) as e:
        ln.examples.cellxgene.create_cellxgene_schema(
            schema_version="5.3.0", field_types=["ensembl_gene_ids"]
        )
    assert "Invalid field_types" in str(e.value)


def test_invalid_schema_Version():
    with pytest.raises(ValueError) as e:
        ln.examples.cellxgene.create_cellxgene_schema(schema_version="200.0.0")
    assert "Invalid schema_version" in str(e.value)
