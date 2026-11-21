import lamindb as ln
import pandas as pd
import pyarrow.parquet as pq
import pytest


def test_parquet_kwargs():
    df = pd.DataFrame(
        {
            "a": [3, 1, 4, 2],
            "b": ["c", "a", "d", "b"],
            "c": [3.3, 1.1, 4.4, 2.2],
        }
    )
    df_sorted = df.sort_values(by=["a", "b"])
    sorting_columns = [
        pq.SortingColumn(0, descending=False, nulls_first=False),
        pq.SortingColumn(1, descending=False, nulls_first=False),
    ]
    artifact = ln.Artifact.from_dataframe(
        df_sorted,
        key="df_sorted.parquet",
        parquet_kwargs={"sorting_columns": sorting_columns},
    ).save()
    pyarrow_dataset = artifact.open()
    fragment = next(pyarrow_dataset.get_fragments())
    assert list(fragment.metadata.row_group(0).sorting_columns) == sorting_columns


def test_curator_refuses_metadata_after_data_changes():
    feature = ln.Feature(name="mutable_sort", dtype=int).save()
    schema = ln.Schema(features=[feature], sorted_by=["mutable_sort"]).save()
    df = pd.DataFrame({"mutable_sort": [1, 2]})
    curator = ln.curators.DataFrameCurator(df, schema)
    curator.validate()
    curator.dataset.iloc[0, 0] = 3
    with pytest.raises(ln.errors.ValidationError, match="not sorted"):
        curator.save_artifact(
            key="mutated-sort.parquet", write_parquet_sorting_metadata=True
        )


def test_metadata_opt_in_requires_sorted_schema():
    with pytest.raises(ln.errors.InvalidArgument, match="schema with sorted_by"):
        ln.Artifact.from_dataframe(
            pd.DataFrame({"a": [1, 2]}),
            key="no-sorted-schema.parquet",
            write_parquet_sorting_metadata=True,
        )


@pytest.mark.parametrize("descending", [False, True])
def test_sorted_schema_metadata_opt_in_and_nulls(descending):
    feature = ln.Feature(name="nullable_sort", dtype=str, nullable=True).save()
    schema = ln.Schema(
        features=[feature], sorted_by=[("nullable_sort", not descending)]
    ).save()
    df = pd.DataFrame(
        {"nullable_sort": ["b", "a", None] if descending else ["a", "b", None]}
    )
    default = ln.Artifact.from_dataframe(
        df, key=f"sort-default-{descending}.parquet", schema=schema
    ).save()
    assert not pq.ParquetFile(default.path).metadata.row_group(0).sorting_columns
    artifact = ln.Artifact.from_dataframe(
        df,
        key=f"sort-opt-in-{descending}.parquet",
        schema=schema,
        write_parquet_sorting_metadata=True,
    ).save()
    assert list(
        pq.ParquetFile(artifact.path).metadata.row_group(0).sorting_columns
    ) == [pq.SortingColumn(0, descending=descending, nulls_first=False)]
    default.delete(permanent=True)
    artifact.delete(permanent=True)
    with pytest.raises(ln.errors.ValidationError, match="not sorted"):
        ln.Artifact.from_dataframe(
            pd.DataFrame({"nullable_sort": [None, "a", "b"]}),
            key=f"sort-invalid-{descending}.parquet",
            schema=schema,
            write_parquet_sorting_metadata=True,
        )
    with pytest.raises(ln.errors.ValidationError, match="not sorted"):
        ln.Artifact.from_dataframe(
            pd.DataFrame({"nullable_sort": [None, "a", "b"]}),
            key=f"sort-invalid-default-{descending}.parquet",
            schema=schema,
        )
