import anndata as ad
import bionty as bt
import lamindb as ln
import numpy as np
import pandas as pd
from lamindb.core import datasets
from lamindb.core._data import _describe_postgres


def check_df_equality(actual_df: pd.DataFrame, expected_df: pd.DataFrame):
    """Checks equality between two DataFrames.

    Special handling for columns containing sets and NaN values.
    """
    pd.testing.assert_index_equal(actual_df.index, expected_df.index)
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
                assert (
                    actual_val == expected_val
                ), f"Set mismatch at index {idx} in column {col}"
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
    ln.Feature(name="cell_medium", dtype="cat[ULabel]").save()
    ln.Feature(name="sample_note", dtype="str").save()
    ln.Feature(name="cell_type_by_expert", dtype="cat[bionty.CellType]").save()
    ln.Feature(name="cell_type_by_model", dtype="cat[bionty.CellType]").save()
    # dataset-level metadata
    ln.Feature(name="temperature", dtype="float").save()
    ln.Feature(name="study", dtype="cat[ULabel]").save()
    ln.Feature(name="date_of_study", dtype="date").save()
    ln.Feature(name="study_note", dtype="str").save()
    ## Permissible values for categoricals
    ln.ULabel.from_values(["DMSO", "IFNG"], create=True).save()
    ln.ULabel.from_values(
        ["Candidate marker study 1", "Candidate marker study 2"], create=True
    ).save()
    bt.CellType.from_values(["B cell", "T cell"], create=True).save()

    ## Ingest dataset1
    adata = datasets.small_dataset1(format="anndata")
    curator = ln.Curator.from_anndata(
        adata,
        var_index=bt.Gene.symbol,
        categoricals={
            "cell_medium": ln.ULabel.name,
            "cell_type_by_expert": bt.CellType.name,
            "cell_type_by_model": bt.CellType.name,
        },
        organism="human",
    )
    artifact = curator.save_artifact(key="example_datasets/dataset1.h5ad")
    artifact.features.add_values(adata.uns)

    # Ingest dataset2
    adata2 = datasets.small_dataset2(format="anndata")
    curator = ln.Curator.from_anndata(
        adata2,
        var_index=bt.Gene.symbol,
        categoricals={
            "cell_medium": ln.ULabel.name,
            "cell_type_by_model": bt.CellType.name,
        },
        organism="human",
    )
    artifact2 = curator.save_artifact(key="example_datasets/dataset2.h5ad")
    artifact2.features.add_values(adata2.uns)

    # Test df(include=[...])
    df = ln.Artifact.df(include=["feature_sets__hash", "feature_sets__name"]).drop(
        "uid", axis=1
    )
    expected_data = {
        "id": [2, 1],
        "key": ["example_datasets/dataset2.h5ad", "example_datasets/dataset1.h5ad"],
        "description": [None, None],
        "feature_sets__hash": [
            set(artifact2.feature_sets.all().values_list("hash", flat=True)),
            set(artifact.feature_sets.all().values_list("hash", flat=True)),
        ],
        "feature_sets__name": [{None}, {None}],
    }
    expected_df = pd.DataFrame(expected_data).set_index("id")
    check_df_equality(df, expected_df)

    # Test df(features=True)
    df = ln.Artifact.df(features=True).drop("uid", axis=1)
    expected_data = {
        "id": [2, 1],
        "key": ["example_datasets/dataset2.h5ad", "example_datasets/dataset1.h5ad"],
        "description": [None, None],
        "cell_type_by_expert": [np.nan, {"T cell", "B cell"}],
        "cell_type_by_model": [{"T cell", "B cell"}, {"T cell", "B cell"}],
        "study": [{"Candidate marker study 2"}, {"Candidate marker study 1"}],
        "cell_medium": [{"IFNG", "DMSO"}, {"IFNG", "DMSO"}],
        "temperature": [{21.6}, np.nan],
        "study_note": [
            {
                "We had a great time performing this study and the results look compelling."
            },
            np.nan,
        ],
        "date_of_study": [{"2024-12-01"}, np.nan],
    }
    expected_df = pd.DataFrame(expected_data).set_index("id")
    check_df_equality(df, expected_df)

    # expected output has italicized elements that can't be tested
    # hence testing is restricted to section content, not headings
    description_tree = _describe_postgres(artifact, print_types=True)

    # general section
    assert len(description_tree.children) == 3
    gernal_node = description_tree.children[0]
    assert gernal_node.label.plain == "General"
    assert gernal_node.children[0].label == f".uid = '{artifact.uid}'"
    assert gernal_node.children[1].label == ".key = 'example_datasets/dataset1.h5ad'"
    assert ".size = " in gernal_node.children[2].label
    assert ".hash = " in gernal_node.children[3].label
    assert gernal_node.children[4].label.plain == ".n_observations = 3"
    assert ".path = " in gernal_node.children[5].label.plain
    assert ".created_by = " in gernal_node.children[6].label.plain
    assert ".created_at = " in gernal_node.children[7].label.plain

    # dataset section
    dataset_node = description_tree.children[1]
    assert dataset_node.label.plain == "Dataset/.feature_sets"
    assert len(dataset_node.children) == 2
    assert len(dataset_node.children[0].label.rows) == 3
    assert len(dataset_node.children[0].label.columns) == 3
    assert dataset_node.children[0].label.columns[0].header.plain == "var • 3"
    assert dataset_node.children[0].label.columns[0]._cells == ["CD8A", "CD4", "CD14"]
    assert dataset_node.children[0].label.columns[1].header.plain == "[bionty.Gene]"
    assert dataset_node.children[0].label.columns[1]._cells[0].plain == "int"
    assert dataset_node.children[1].label.columns[0].header.plain == "obs • 4"
    assert dataset_node.children[1].label.columns[0]._cells == [
        "cell_medium",
        "cell_type_by_expert",
        "cell_type_by_model",
    ]
    assert dataset_node.children[1].label.columns[1].header.plain == "[Feature]"
    assert dataset_node.children[1].label.columns[1]._cells[0].plain == "cat[ULabel]"
    assert (
        dataset_node.children[1].label.columns[1]._cells[1].plain
        == "cat[bionty.CellType]"
    )
    assert (
        dataset_node.children[1].label.columns[1]._cells[2].plain
        == "cat[bionty.CellType]"
    )
    assert dataset_node.children[1].label.columns[2]._cells == [
        "DMSO, IFNG",
        "B cell, T cell",
        "B cell, T cell",
    ]

    # annotations section
    annotations_node = description_tree.children[2]
    assert annotations_node.label.plain == "Annotations"
    assert len(annotations_node.children) == 2
    assert len(annotations_node.children[0].label.columns) == 3
    assert len(annotations_node.children[0].label.rows) == 4
    assert annotations_node.children[0].label.columns[0].header.plain == "Features"
    assert annotations_node.children[0].label.columns[0]._cells == [
        "study",
        "date_of_study",
        "study_note",
        "temperature",
    ]
    assert (
        annotations_node.children[0].label.columns[1]._cells[0].plain == "cat[ULabel]"
    )
    assert annotations_node.children[0].label.columns[1]._cells[1].plain == "date"
    assert annotations_node.children[0].label.columns[1]._cells[2].plain == "str"
    assert annotations_node.children[0].label.columns[1]._cells[3].plain == "float"
    assert annotations_node.children[0].label.columns[2]._cells == [
        "Candidate marker study 1",
        "2024-12-01",
        "We had a great time performing this study and the results look compelling.",
        "21.6",
    ]
    assert len(annotations_node.children[0].label.columns) == 3
    assert len(annotations_node.children[1].label.rows) == 2
    assert annotations_node.children[1].label.columns[0].header.plain == "Labels"
    assert annotations_node.children[1].label.columns[0]._cells == [
        ".cell_types",
        ".ulabels",
    ]
    assert (
        annotations_node.children[1].label.columns[1]._cells[0].plain
        == "bionty.CellType"
    )
    assert annotations_node.children[1].label.columns[1]._cells[1].plain == "ULabel"
    assert annotations_node.children[1].label.columns[2]._cells == [
        "'B cell', 'T cell'",
        "'DMSO', 'IFNG', 'Candidate marker study 1'",
    ]

    artifact.delete(permanent=True)
    artifact2.delete(permanent=True)
    ln.FeatureSet.filter().delete()
    ln.Feature.filter().delete()
    bt.Gene.filter().delete()
    ln.ULabel.filter().delete()
    bt.CellType.filter().delete()
