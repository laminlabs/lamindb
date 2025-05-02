from __future__ import annotations

from typing import Any, Literal

import anndata as ad
import numpy as np
import pandas as pd


def small_dataset3_cellxgene(
    otype: Literal["DataFrame", "AnnData"] = "AnnData",
) -> tuple[pd.DataFrame, dict[str, Any]] | ad.AnnData:
    # TODO: consider other ids for other organisms
    # "ENSMUSG00002076988"
    var_ids = ["invalid_ensembl_id", "ENSG00000000419", "ENSG00000139618"]
    dataset_dict = {
        var_ids[0]: [2, 3, 3],
        var_ids[1]: [3, 4, 5],
        var_ids[2]: [4, 2, 3],
        "disease_ontology_term_id": ["MONDO:0004975", "MONDO:0004980", "MONDO:0004980"],
        "organism": ["human", "human", "human"],
        "sex": ["female", "male", "unknown"],
        "sex_ontology_term_id": ["PATO:0000383", "PATO:0000384", "unknown"],
        "tissue": ["lungg", "lungg", "heart"],
        "donor": ["-1", "1", "2"],
    }
    dataset_df = pd.DataFrame(
        dataset_dict,
        index=["barcode1", "barcode2", "barcode3"],
    )
    dataset_df["tissue"] = dataset_df["tissue"].astype("category")
    ad.AnnData(
        dataset_df[var_ids],
        obs=dataset_df[[key for key in dataset_dict if key not in var_ids]],
    )
    if otype == "DataFrame":
        return dataset_df
    else:
        dataset_ad = ad.AnnData(dataset_df.iloc[:, :3], obs=dataset_df.iloc[:, 3:])
        return dataset_ad


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
