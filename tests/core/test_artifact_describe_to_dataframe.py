import bionty as bt
import lamindb as ln
import numpy as np
import pandas as pd
import pytest
from lamindb.models._describe import describe_postgres


def _check_df_equality(actual_df: pd.DataFrame, expected_df: pd.DataFrame) -> bool:
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


def test_describe_artifact_from_remote_instance(capsys):
    # test describing from a remote instance with less modules
    artifact = ln.Artifact.connect("laminlabs/lamin-site-assets").first()
    artifact.describe()
    captured = capsys.readouterr()
    assert len(captured.out) > 50
    assert "artifact" in captured.out.lower()


# parallels the `registries` guide
# please also see the test_querset.py tests
def test_describe_to_dataframe_example_dataset():
    ln.examples.datasets.mini_immuno.save_mini_immuno_datasets()
    artifact = ln.Artifact.get(key="examples/dataset1.h5ad")
    artifact2 = ln.Artifact.get(key="examples/dataset2.h5ad")

    with pytest.raises(ValueError) as error:
        artifact.features.remove_values("cell_type_by_expert")
    assert "Cannot remove values for dataset features." in error.exconly()

    # Test df(include=[...])
    df = (
        ln.Artifact.filter(key__startswith="examples/dataset", suffix=".h5ad")
        .order_by("-key")
        .to_dataframe(include=["feature_sets__hash", "feature_sets__name"])
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
    _check_df_equality(df, expected_df)

    # Test df with features
    # test that the records filter DOES NOT affect joining the annotations
    # we want it to only affect the artifact query (even though here, it won't change the result as both artifacts have the IFNG label)
    df = (
        ln.Artifact.filter(
            key__startswith="examples/dataset",
            suffix=".h5ad",
            records__name="IFNG",
        )
        .order_by("-key")
        .to_dataframe(
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
    _check_df_equality(df, expected_df)

    # expected output has italicized elements that can't be tested
    # hence testing is restricted to section content, not headings
    output = artifact.describe(return_str=True)
    assert "hash:" in output
    assert "size:" in output
    assert "n_observations: 3" in output
    assert "storage/path:" in output
    assert "created_by:" in output
    assert "created_at:" in output

    # dataset section
    assert (
        artifact.features.describe(return_str=True)
        == """Artifact: examples/dataset1.h5ad (0000)
├── Dataset features
│   ├── obs (4)
│   │   cell_type_by_expe…  bionty.CellType         B cell, CD8-positive, alpha…
│   │   cell_type_by_model  bionty.CellType         B cell, T cell
│   │   perturbation        Record                  DMSO, IFNG
│   │   sample_note         str
│   └── var.T (3 bionty.G…
│       CD8A                num
│       CD4                 num
│       CD14                num
└── External features
    └── experiment          Record                  Experiment 1
        date_of_study       date                    2024-12-01
        study_metadata      dict                    {'detail1': '123', 'detail2…
        study_note          str                     We had a great time perform…
        temperature         float                   21.6"""
    )

    # labels section
    description_tree = describe_postgres(artifact)
    labels_node = description_tree.children[-1].label
    assert labels_node.label.plain == "Labels"
    assert len(labels_node.children[0].label.columns) == 3
    assert len(labels_node.children[0].label.rows) == 2
    assert labels_node.children[0].label.columns[0]._cells == [
        ".records",
        ".cell_types",
    ]
    assert labels_node.children[0].label.columns[1]._cells[0].plain == "Record"
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

    # test that only external feature are removed upon artifact.features.remove_values()
    before = artifact.features.get_values()
    adata = artifact.load()
    just_internal = {}
    for col in adata.obs.columns:
        if col in before:
            just_internal[col] = before[col]
    artifact.features.remove_values()
    assert just_internal == artifact.features.get_values()

    artifact.delete(permanent=True)
    artifact2.delete(permanent=True)
    ln.Schema.get(name="anndata_ensembl_gene_ids_and_valid_features_in_obs").delete(
        permanent=True
    )
    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    bt.Gene.filter().delete(permanent=True)
    ln.Record.filter().delete(permanent=True)
    bt.CellType.filter().delete(permanent=True)
