import bionty as bt
import lamindb as ln
import numpy as np
import pandas as pd
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
    ln.examples.ingest_mini_immuno_datasets()
    artifact = ln.Artifact.get(key="examples/dataset1.h5ad")
    artifact2 = ln.Artifact.get(key="examples/dataset2.h5ad")

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
        "experiment": pd.Categorical(["Experiment 2", "Experiment 1"]),
        "perturbation": [{"IFNG", "DMSO"}, {"IFNG", "DMSO"}],
        "temperature": [22.6, 21.6],
        "study_note": [
            np.nan,
            "We had a great time performing this study and the results look compelling.",
        ],
        "date_of_study": ["2025-02-13", "2024-12-01"],
        "study_metadata": [
            {"detail1": "456", "detail2": 2},
            {"detail1": "123", "detail2": 1},
        ],
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

    # The structure is now different due to two-column layout
    # First few children are two-column pairs, then single-column items

    # Check that uid appears in the first two-column row
    first_row = general_node.children[0].label.plain
    assert f"key: {artifact.key}" in first_row

    # Check that hash appears somewhere in the two-column section
    found_hash = False
    found_size = False
    found_n_observations = False

    # Look through the two-column rows for hash, size, and n_observations
    for child in general_node.children:
        child_text = child.label.plain
        if "hash: " in child_text:
            found_hash = True
        if "size: " in child_text:
            found_size = True
        if "n_observations: 3" in child_text:
            found_n_observations = True

    assert found_hash, "Hash should be present in the general section"
    assert found_size, "Size should be present in the general section"
    assert found_n_observations, (
        "n_observations should be present in the general section"
    )

    # Check single-column items (these appear after the two-column items)
    found_key = False
    found_path = False
    found_created_by = False
    found_created_at = False

    for child in general_node.children:
        child_text = child.label.plain
        if "key: examples/dataset1.h5ad" in child_text:
            found_key = True
        if "storage path: " in child_text:
            found_path = True
        if "created_by: " in child_text:
            found_created_by = True
        if "created_at: " in child_text:
            found_created_at = True

    assert found_key, "Key should be present in the general section"
    assert found_path, "Storage path should be present in the general section"
    assert found_created_by, "Created by should be present in the general section"
    assert found_created_at, "Created at should be present in the general section"

    # dataset section
    assert (
        artifact.features.describe(return_str=True)
        == """Artifact .h5ad · AnnData · dataset
├── Dataset features
│   ├── obs • 4             [Feature]
│   │   cell_type_by_expe…  cat[bionty.CellType]    B cell, CD8-positive, alpha…
│   │   cell_type_by_model  cat[bionty.CellType]    B cell, T cell
│   │   perturbation        cat[ULabel]             DMSO, IFNG
│   │   sample_note         str
│   └── var.T • 3           [bionty.Gene.ensembl_…
│       CD8A                num
│       CD4                 num
│       CD14                num
└── Linked features
    └── experiment          cat[ULabel]             Experiment 1
        date_of_study       date                    2024-12-01
        study_metadata      dict                    {'detail1': '123', 'detail2…
        study_note          str                     We had a great time perform…
        temperature         float                   21.6"""
    )

    # labels section
    labels_node = description_tree.children[3].label
    assert labels_node.label.plain == "Labels"
    assert len(labels_node.children[0].label.columns) == 3
    assert len(labels_node.children[0].label.rows) == 2
    assert labels_node.children[0].label.columns[0]._cells == [
        ".ulabels",
        ".cell_types",
    ]
    assert labels_node.children[0].label.columns[1]._cells[0].plain == "ULabel"
    assert labels_node.children[0].label.columns[1]._cells[1].plain == "bionty.CellType"
    assert {
        c.strip()
        for c in ",".join(labels_node.children[0].label.columns[2]._cells).split(",")
    } == {
        "DMSO",
        "IFNG",
        "Experiment 1",
        "B cell",
        "T cell",
        "CD8-positive",
        "alpha-beta T cell",
    }

    artifact.delete(permanent=True)
    artifact2.delete(permanent=True)
    ln.Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs").delete()
    ln.Schema.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()
