import lamindb as ln
import pandas as pd
import pyarrow.parquet as pq


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
