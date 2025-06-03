"""The mini immuno dataset.

.. autosummary::
   :toctree: .

   define_features_labels
   get_dataset1
   get_dataset2

"""

from __future__ import annotations

from typing import TYPE_CHECKING, Literal

import anndata as ad
import pandas as pd

if TYPE_CHECKING:
    from lamindb.models import Schema


def define_features_labels() -> None:
    """Features & labels to validate the mini immuno datasets.

    .. literalinclude:: scripts/define_mini_immuno_features_labels.py
        :language: python
    """
    import sys
    from pathlib import Path

    docs_path = Path(__file__).parent.parent.parent.parent / "docs" / "scripts"
    if str(docs_path) not in sys.path:
        sys.path.append(str(docs_path))

    import define_mini_immuno_features_labels  # noqa


def define_mini_immuno_schema_flexible() -> Schema:
    """Features & labels to validate the mini immuno datasets.

    .. literalinclude:: scripts/define_mini_immuno_schema_flexible.py
        :language: python
    """
    import sys
    from pathlib import Path

    from lamindb.models import Schema

    docs_path = Path(__file__).parent.parent.parent.parent / "docs" / "scripts"
    if str(docs_path) not in sys.path:
        sys.path.append(str(docs_path))

    define_features_labels()
    import define_mini_immuno_schema_flexible  # noqa

    return Schema.get(name="Mini immuno schema")


def get_dataset1(
    otype: Literal["DataFrame", "AnnData"] = "DataFrame",
    gene_symbols_in_index: bool = False,
    with_typo: bool = False,
    with_cell_type_synonym: bool = False,
    with_cell_type_typo: bool = False,
    with_gene_typo: bool = False,
    with_outdated_gene: bool = False,
    with_wrong_subtype: bool = False,
    with_index_type_mismatch: bool = False,
) -> pd.DataFrame | ad.AnnData:
    """A small tabular dataset measuring expression & metadata."""
    # define the data in the dataset
    # it's a mix of numerical measurements and observation-level metadata
    ifng = "IFNJ" if with_typo else "IFNG"
    thing = "ulabel_but_not_perturbation" if with_wrong_subtype else "DMSO"
    if gene_symbols_in_index:
        var_ids = ["CD8A", "CD4", "CD14" if not with_gene_typo else "GeneTypo"]
    else:
        var_ids = [
            "ENSG00000153563",
            "ENSG00000010610",
            "ENSG00000170458"
            if not with_gene_typo
            else "GeneTypo"
            if not with_outdated_gene
            else "ENSG00000278198",
        ]
    abt_cell = (
        "CD8-pos alpha-beta T cell"
        if with_cell_type_typo
        else "CD8-positive, alpha-beta T cell"
    )
    dataset_dict = {
        var_ids[0]: [1, 2, 3],
        var_ids[1]: [3, 4, 5],
        var_ids[2]: [5, 6, 7],
        "perturbation": pd.Categorical(["DMSO", ifng, thing]),
        "sample_note": ["was ok", "looks naah", "pretty! 🤩"],
        "cell_type_by_expert": pd.Categorical(
            ["B-cell" if with_cell_type_synonym else "B cell", abt_cell, abt_cell]
        ),
        "cell_type_by_model": pd.Categorical(["B cell", "T cell", "T cell"]),
        "assay_oid": pd.Categorical(["EFO:0008913", "EFO:0008913", "EFO:0008913"]),
        "concentration": ["0.1%", "200 nM", "0.1%"],
        "treatment_time_h": [24, 24, 6],
        "donor": ["D0001", "D0002", None],
        "donor_ethnicity": [
            ["Chinese", "Singaporean Chinese"],
            ["Chinese", "Han Chinese"],
            ["Chinese"],
        ],
    }
    # define the dataset-level metadata
    metadata = {
        "temperature": 21.6,
        "experiment": "Experiment 1",
        "date_of_study": "2024-12-01",
        "study_note": "We had a great time performing this study and the results look compelling.",
    }
    # the dataset as DataFrame
    dataset_df = pd.DataFrame(
        dataset_dict,
        index=["sample1", "sample2", 0]  # type: ignore
        if with_index_type_mismatch
        else ["sample1", "sample2", "sample3"],
    )
    if otype == "DataFrame":
        for key, value in metadata.items():
            dataset_df.attrs[key] = value
        return dataset_df
    else:
        del dataset_df[
            "donor_ethnicity"
        ]  # remove the donor_ethnicity because AnnData save will error
        dataset_ad = ad.AnnData(
            dataset_df.iloc[:, :3], obs=dataset_df.iloc[:, 3:], uns=metadata
        )
        return dataset_ad


def get_dataset2(
    otype: Literal["DataFrame", "AnnData"],
    gene_symbols_in_index: bool = False,
) -> pd.DataFrame | ad.AnnData:
    if gene_symbols_in_index:
        var_ids = ["CD8A", "CD4", "CD38"]
    else:
        var_ids = ["ENSG00000153563", "ENSG00000010610", "ENSG00000004468"]
    dataset_dict = {
        var_ids[0]: [2, 3, 3],
        var_ids[1]: [3, 4, 5],
        var_ids[2]: [4, 2, 3],
        "perturbation": pd.Categorical(["DMSO", "IFNG", "IFNG"]),
        "cell_type_by_model": pd.Categorical(["B cell", "T cell", "T cell"]),
        "concentration": ["0.1%", "200 nM", "0.1%"],
        "treatment_time_h": [24, 24, 6],
        "donor": ["D0003", "D0003", "D0004"],
    }
    metadata = {
        "temperature": 22.6,
        "experiment": "Experiment 2",
        "date_of_study": "2025-02-13",
    }
    dataset_df = pd.DataFrame(
        dataset_dict,
        index=["sample4", "sample5", "sample6"],
    )
    ad.AnnData(
        dataset_df[var_ids],
        obs=dataset_df[["perturbation", "cell_type_by_model"]],
    )
    if otype == "DataFrame":
        for key, value in metadata.items():
            dataset_df.attrs[key] = value
        return dataset_df
    else:
        dataset_ad = ad.AnnData(
            dataset_df.iloc[:, :3], obs=dataset_df.iloc[:, 3:], uns=metadata
        )
        return dataset_ad
