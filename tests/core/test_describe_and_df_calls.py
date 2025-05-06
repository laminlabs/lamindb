import bionty as bt
import lamindb as ln
import numpy as np
import pandas as pd
from lamindb.core import datasets
from lamindb.models.artifact import _describe_postgres


def check_df_equality(actual_df: pd.DataFrame, expected_df: pd.DataFrame):
    """Checks equality between two DataFrames.

    Special handling for columns containing sets and NaN values.
    """
    # do not test indices by default
    # pd.testing.assert_index_equal(actual_df.index, expected_df.index)
    expected_df.index = actual_df.index
    assert set(actual_df.columns) == set(expected_df.columns)
    for col in expected_df.columns:
        # Detect if column contains sets by checking first non-null value
        first_value = next((v for v in expected_df[col] if pd.notna(v)), None)
        is_set_column = isinstance(first_value, set)
        if is_set_column:
            # For set columns, compare sets with NaN handling
            for idx in expected_df.index:
                actual_val = actual_df.loc[idx, col]
                expected_val = expected_df.loc[idx, col]
                # If both are NaN, they're equal
                if pd.isna(actual_val) and pd.isna(expected_val):
                    continue
                # If one is NaN and the other isn't, they're not equal
                if pd.isna(actual_val) != pd.isna(expected_val):
                    raise AssertionError(f"NaN mismatch at index {idx} in column {col}")
                # If neither is NaN, compare the sets
                assert actual_val == expected_val, (
                    f"Set mismatch at index {idx} in column {col}"
                )
        else:
            pd.testing.assert_series_equal(
                actual_df[col],
                expected_df[col],
                check_names=False,  # ignore series names
            )
    return True


# parallels the `registries` guide
# please also see the test_querset.py tests
def test_curate_df():
    ## Create a more complex case
    # observation-level metadata
    ln.Feature(name="perturbation", dtype="cat[ULabel]").save()
    ln.Feature(name="sample_note", dtype="str").save()
    ln.Feature(name="cell_type_by_expert", dtype="cat[bionty.CellType]").save()
    ln.Feature(name="cell_type_by_model", dtype="cat[bionty.CellType]").save()
    # dataset-level metadata
    ln.Feature(name="temperature", dtype="float").save()
    ln.Feature(name="experiment", dtype="cat[ULabel]").save()
    ln.Feature(name="date_of_study", dtype="date").save()
    ln.Feature(name="study_note", dtype="str").save()
    ln.Feature(name="study_metadata", dtype=dict).save()

    ## Permissible values for categoricals
    ln.ULabel.from_values(["DMSO", "IFNG"], create=True).save()
    ln.ULabel.from_values(["Experiment 1", "Experiment 2"], create=True).save()
    bt.CellType.from_values(["B cell", "T cell"], create=True).save()

    ## Ingest dataset1
    adata = datasets.small_dataset1(otype="AnnData")
    curator = ln.Curator.from_anndata(
        adata,
        var_index=bt.Gene.ensembl_gene_id,
        categoricals={
            "perturbation": ln.ULabel.name,
            "cell_type_by_expert": bt.CellType.name,
            "cell_type_by_model": bt.CellType.name,
        },
    )
    artifact = curator.save_artifact(key="examples/dataset1.h5ad")
    d1 = {"study_metadata": {"a": "123", "b": 1}}
    dataset_metadata = adata.uns
    dataset_metadata.update(d1)
    artifact.features.add_values(dataset_metadata)
    print("dataset1", adata.uns)

    # Ingest dataset2
    adata2 = datasets.small_dataset2(otype="AnnData")
    curator = ln.Curator.from_anndata(
        adata2,
        var_index=bt.Gene.ensembl_gene_id,
        categoricals={
            "perturbation": ln.ULabel.name,
            "cell_type_by_model": bt.CellType.name,
        },
    )
    artifact2 = curator.save_artifact(key="examples/dataset2.h5ad")
    d2 = {"study_metadata": {"a": "456", "b": 2}}
    dataset_metadata = adata2.uns
    dataset_metadata.update(d2)
    artifact2.features.add_values(dataset_metadata)
    print("dataset2", adata2.uns)

    # Test df(include=[...])
    df = (
        ln.Artifact.filter(key__startswith="examples/dataset", suffix=".h5ad")
        .order_by("-key")
        .df(include=["feature_sets__hash", "feature_sets__name"])
        .drop(["uid"], axis=1)
    )
    expected_data = {
        "key": ["examples/dataset2.h5ad", "examples/dataset1.h5ad"],
        "feature_sets__hash": [
            set(artifact2.feature_sets.all().values_list("hash", flat=True)),
            set(artifact.feature_sets.all().values_list("hash", flat=True)),
        ],
        "feature_sets__name": [{None}, {None}],
    }
    expected_df = pd.DataFrame(expected_data)
    check_df_equality(df, expected_df)

    # Test df with features
    # test that the ulabels filter DOES NOT affect joining the annotations
    # we want it to only affect the artifact query (even though here, it won't change the result as both artifacts have the IFNG label)
    df = (
        ln.Artifact.filter(
            key__startswith="examples/dataset",
            suffix=".h5ad",
            ulabels__name="IFNG",
        )
        .order_by("-key")
        .df(
            features=[
                "cell_type_by_expert",
                "cell_type_by_model",
                "experiment",
                "perturbation",
                "temperature",
                "study_note",
                "date_of_study",
            ]
        )
        .drop(["uid"], axis=1)
    )
    expected_data = {
        "key": ["examples/dataset2.h5ad", "examples/dataset1.h5ad"],
        "cell_type_by_expert": [np.nan, {"CD8-positive, alpha-beta T cell", "B cell"}],
        "cell_type_by_model": [{"T cell", "B cell"}, {"T cell", "B cell"}],
        "experiment": ["Experiment 2", "Experiment 1"],
        "perturbation": [{"IFNG", "DMSO"}, {"IFNG", "DMSO"}],
        "temperature": [22.6, 21.6],
        "study_note": [
            np.nan,
            "We had a great time performing this study and the results look compelling.",
        ],
        "date_of_study": ["2025-02-13", "2024-12-01"],
        "study_metadata": [{"a": "456", "b": 2}, {"a": "123", "b": 1}],
    }
    expected_df = pd.DataFrame(expected_data)
    check_df_equality(df, expected_df)

    # expected output has italicized elements that can't be tested
    # hence testing is restricted to section content, not headings
    description_tree = _describe_postgres(artifact)

    # general section
    assert (
        len(description_tree.children) == 4
    )  # general, internal features, external features, labels
    general_node = description_tree.children[0]
    assert general_node.label.plain == "General"
    assert general_node.children[0].label == f".uid = '{artifact.uid}'"
    assert general_node.children[1].label == ".key = 'examples/dataset1.h5ad'"
    assert ".size = " in general_node.children[2].label
    assert ".hash = " in general_node.children[3].label
    assert general_node.children[4].label.plain == ".n_observations = 3"
    assert ".path = " in general_node.children[5].label.plain
    assert ".created_by = " in general_node.children[6].label.plain
    assert ".created_at = " in general_node.children[7].label.plain

    # dataset section
    assert (
        artifact.features.describe(return_str=True)
        == """Artifact .h5ad/AnnData
├── Dataset features
│   ├── var • 3             [bionty.Gene]
│   │   CD8A                int
│   │   CD4                 int
│   │   CD14                int
│   └── obs • 4             [Feature]
│       cell_type_by_expe…  cat[bionty.CellT…  B cell, CD8-positive, alpha-beta…
│       cell_type_by_model  cat[bionty.CellT…  B cell, T cell
│       perturbation        cat[ULabel]        DMSO, IFNG
│       sample_note         str
└── Linked features
    └── experiment          cat[ULabel]        Experiment 1
        date_of_study       date               2024-12-01
        study_metadata      dict               {'a': '123', 'b': 1}
        study_note          str                We had a great time performing t…
        temperature         float              21.6"""
    )

    # labels section
    labels_node = description_tree.children[3].label
    assert labels_node.label.plain == "Labels"
    assert len(labels_node.children[0].label.columns) == 3
    assert len(labels_node.children[0].label.rows) == 2
    assert labels_node.children[0].label.columns[0]._cells == [
        ".cell_types",
        ".ulabels",
    ]
    assert labels_node.children[0].label.columns[1]._cells[0].plain == "bionty.CellType"
    assert labels_node.children[0].label.columns[1]._cells[1].plain == "ULabel"
    assert labels_node.children[0].label.columns[2]._cells == [
        "B cell, T cell, CD8-positive, alpha-beta T cell",
        "DMSO, IFNG, Experiment 1",
    ]

    artifact.delete(permanent=True)
    artifact2.delete(permanent=True)
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()
