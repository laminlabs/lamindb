import lamindb as ln


def test_cxg_curator():
    schema_version = "5.2.0"
    organism = "human"
    adata = ln.core.datasets.small_dataset3_cellxgene()
    curator = ln.curators.CellxGeneAnnDataCatManager(
        adata, organism=organism, schema_version=schema_version
    )
    adata.obs.rename(columns={"donor": "donor_id"}, inplace=True)
    curator = ln.curators.CellxGeneAnnDataCatManager(
        adata,
        defaults=ln.curators.CellxGeneAnnDataCatManager._get_categoricals_defaults(),
        organism=organism,
        schema_version=schema_version,
    )
    assert not curator.validate()
    adata = adata[:, ~adata.var.index.isin(curator.non_validated["var_index"])]
    adata.obs["tissue"] = adata.obs["tissue"].cat.rename_categories({"lungg": "lung"})
    curator = ln.curators.CellxGeneAnnDataCatManager(
        adata, organism=organism, schema_version=schema_version
    )
    assert curator.validate()
