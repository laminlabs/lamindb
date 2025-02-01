import lamindb as ln


def test_cxg_curator():
    adata = ln.core.datasets.small_dataset3_cellxgene()
    curator = ln.curators.CellxGeneAnnDataCatCurator(
        adata, organism="human", schema_version="5.1.0"
    )
    adata.obs.rename(columns={"donor": "donor_id"}, inplace=True)
    curator = ln.curators.CellxGeneAnnDataCatCurator(
        adata,
        defaults=ln.curators.CellxGeneAnnDataCatCurator._get_categoricals_defaults(),
        organism="human",
        schema_version="5.1.0",
    )
    assert not curator.validate()
    adata = adata[:, ~adata.var.index.isin(curator.non_validated["var_index"])]
    adata.obs["tissue"] = adata.obs["tissue"].cat.rename_categories({"lungg": "lung"})
    curator = ln.curators.CellxGeneAnnDataCatCurator(
        adata, organism="human", schema_version="5.1.0"
    )
    assert curator.validate()
