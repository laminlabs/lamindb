from __future__ import annotations

from typing import Any, Literal

import anndata as ad
import numpy as np
import pandas as pd


def small_dataset3_cellxgene(
    otype: Literal["DataFrame", "AnnData"] = "AnnData",
    with_obs_defaults: bool = False,
    with_obs_typo: bool = False,
) -> tuple[pd.DataFrame, dict[str, Any]] | ad.AnnData:
    # TODO: consider other ids for other organisms
    # "ENSMUSG00002076988"
    var_ids = ["invalid_ensembl_id", "ENSG00000000419", "ENSG00000139618"]

    lung_id = "UBERON:0002048XXX" if with_obs_typo else "UBERON:0002048"
    obs_df = pd.DataFrame(
        {
            "disease_ontology_term_id": [
                "MONDO:0004975",
                "MONDO:0004980",
                "MONDO:0004980",
            ],
            "development_stage_ontology_term_id": ["unknown", "unknown", "unknown"],
            "organism": ["human", "human", "human"],
            "sex_ontology_term_id": ["PATO:0000383", "PATO:0000384", "unknown"],
            "tissue_ontology_term_id": [lung_id, lung_id, "UBERON:0000948"],
            "cell_type": ["T cell", "B cell", "B cell"],
            "self_reported_ethnicity": ["South Asian", "South Asian", "South Asian"],
            "donor_id": ["-1", "1", "2"],
            "is_primary_data": [False, False, False],
            "suspension_type": ["cell", "cell", "cell"],
            "tissue_type": ["tissue", "tissue", "tissue"],
        },
        index=["barcode1", "barcode2", "barcode3"],
    )

    var_df = pd.DataFrame(
        index=var_ids, data={"feature_is_filtered": [False, False, False]}
    )

    X = pd.DataFrame(
        {
            var_ids[0]: [2, 3, 3],
            var_ids[1]: [3, 4, 5],
            var_ids[2]: [4, 2, 3],
        },
        index=["barcode1", "barcode2", "barcode3"],
        dtype="float32",
    )

    obs_df["donor_id"] = obs_df["donor_id"].astype("category")

    if otype == "DataFrame":
        return pd.concat([X, obs_df], axis=1)
    else:
        adata = ad.AnnData(X=X, obs=obs_df, var=var_df)
        adata.uns["title"] = "CELLxGENE example"
        adata.obsm["X_pca"] = np.array(
            [[-1.2, 0.8], [0.5, -0.3], [0.7, -0.5]], dtype="float32"
        )
        # CELLxGENE requires the `.raw` slot to be set - https://github.com/chanzuckerberg/single-cell-curation/issues/1304
        adata.raw = adata.copy()
        adata.raw.var.drop(columns="feature_is_filtered", inplace=True)
        if with_obs_defaults:
            adata.obs["assay"] = "single-cell RNA sequencing"
        return adata


def anndata_with_obs() -> ad.AnnData:
    """Create a mini anndata with cell_type, disease and tissue."""
    import anndata as ad
    import bionty.base as bionty_base

    celltypes = ["T cell", "hematopoietic stem cell", "hepatocyte", "my new cell type"]
    celltype_ids = ["CL:0000084", "CL:0000037", "CL:0000182", ""]
    diseases = [
        "chronic kidney disease",
        "liver lymphoma",
        "cardiac ventricle disorder",
        "Alzheimer disease",
    ]
    tissues = ["kidney", "liver", "heart", "brain"]
    df = pd.DataFrame()
    df["cell_type"] = celltypes * 10
    df["cell_type_id"] = celltype_ids * 10
    df["tissue"] = tissues * 10
    df["disease"] = diseases * 10
    df.index = "obs" + df.index.astype(str)

    adata = ad.AnnData(X=np.zeros(shape=(40, 100), dtype=np.float32), obs=df)
    adata.var.index = bionty_base.Gene().df().head(100)["ensembl_gene_id"].values

    return adata
