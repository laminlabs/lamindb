import lamindb as ln
import numpy as np


def test_cxg_curator():
    schema_version = "5.2.0"
    adata = ln.core.datasets.small_dataset3_cellxgene()
    curator = ln.curators._legacy.CellxGeneAnnDataCatManager(
        adata, schema_version=schema_version
    )

    adata.obs.rename(columns={"donor": "donor_id"}, inplace=True)
    curator = ln.curators._legacy.CellxGeneAnnDataCatManager(
        adata,
        defaults=ln.curators._legacy.CellxGeneAnnDataCatManager.cxg_categoricals_defaults,
        schema_version=schema_version,
    )
    assert not curator.validate()

    adata = adata[:, ~adata.var.index.isin(curator.non_validated["var_index"])]
    adata.obs["tissue"] = adata.obs["tissue"].cat.rename_categories({"lungg": "lung"})
    curator = ln.curators._legacy.CellxGeneAnnDataCatManager(
        adata, schema_version=schema_version
    )
    assert curator.validate()

    artifact = curator.save_artifact(
        key=f"examples/dataset-curated-against-cxg-{curator.schema_version}.h5ad"
    )
    title = "Cross-tissue immune cell analysis reveals tissue-specific features in humans (for test demo only)"

    adata.obsm["X_umap"] = np.zeros((adata.shape[0], 2))
    adata_cxg = curator.to_cellxgene_anndata(is_primary_data=True, title=title)
    assert "cell_type_ontology_term_id" in adata_cxg.obs.columns

    artifact.delete(permanent=True)
