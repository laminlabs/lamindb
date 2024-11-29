from __future__ import annotations

from typing import Any, Literal

import anndata as ad
import numpy as np
import pandas as pd


def small_dataset1(
    format: Literal["df", "anndata"],
) -> tuple[pd.DataFrame, dict[str, Any]] | ad.AnnData:
    # define the data in the dataset
    # it's a mix of numerical measurements and observation-level metadata
    dataset_dict = {
        "CD8A": [1, 2, 3],
        "CD4": [3, 4, 5],
        "CD14": [5, 6, 7],
        "cell_medium": ["DMSO", "IFNG", "DMSO"],
        "sample_note": ["was ok", "looks naah", "pretty! 🤩"],
        "cell_type_by_expert": ["B cell", "T cell", "T cell"],
        "cell_type_by_model": ["B cell", "T cell", "T cell"],
    }
    # define the dataset-level metadata
    metadata = {
        "temperature": 21.6,
        "study": "Candidate marker study 1",
        "date_of_study": "2024-12-01",
        "study_note": "We had a great time performing this study and the results look compelling.",
    }
    # the dataset as DataFrame
    dataset_df = pd.DataFrame(dataset_dict, index=["sample1", "sample2", "sample3"])
    dataset_ad = ad.AnnData(
        dataset_df.iloc[:, :3], obs=dataset_df.iloc[:, 3:], uns=metadata
    )
    if format == "df":
        return dataset_df, metadata
    else:
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
