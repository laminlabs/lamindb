"""LanceDB table artifacts."""

import shutil
from pathlib import Path
from typing import cast

import lamindb as ln
import pandas as pd
import pytest

lancedb = pytest.importorskip("lancedb")


def _cells(names: list[str | None]) -> pd.DataFrame:
    return pd.DataFrame(
        {
            "cell_type": names,
            "n_genes": list(range(len(names))),
        }
    )


def test_lance_table_roundtrip_filter_schema_and_versions(tmp_path: Path):
    df = _cells(["T cell", None, "B cell"])
    artifact = ln.Artifact.from_dataframe(df, key="lance/cells.lance").save()
    other = None
    revised = None
    collection = None
    schema = None
    feature = None
    try:
        assert artifact.suffix == ".lance"
        assert artifact.key == "lance/cells.lance"
        assert artifact.otype == "DataFrame"
        assert artifact.kind == "dataset"
        assert artifact.n_observations == 3
        assert artifact.overwrite_versions is False
        assert artifact.path.is_dir()

        loaded = cast(pd.DataFrame, artifact.load())
        assert list(loaded["n_genes"]) == [0, 1, 2]
        with pytest.raises(ValueError, match="open_lance"):
            artifact.open()
        table = artifact.open_lance()
        filtered = (
            table.search()
            .where("cell_type IS NOT NULL", prefilter=True)
            .limit(5)
            .to_pandas()
        )
        assert set(filtered["cell_type"]) == {"T cell", "B cell"}

        copied = tmp_path / "copied.lance"
        # A path that is not yet the saved artifact still reports its row count.
        shutil.copytree(artifact.path, copied)
        from_path = ln.Artifact.from_dataframe(copied, key="lance/copied.lance")
        assert from_path.suffix == ".lance"
        assert from_path.n_observations == 3

        feature = ln.Feature(name="lance_schema_feature", dtype=str).save()
        schema = ln.Schema(features=[feature], suffix=".lance").save()
        assert schema.suffix == ".lance"

        previous_path = artifact.path
        revised = ln.Artifact.from_dataframe(
            _cells(["NK cell", "T cell"]), key="lance/cells.lance"
        ).save()
        assert previous_path.exists()
        assert revised.path != previous_path
        assert revised.n_observations == 2
        assert revised.overwrite_versions is False
        assert artifact.open_lance().count_rows() == 3

        other = ln.Artifact.from_dataframe(
            _cells(["B cell"]), key="lance/other.lance"
        ).save()
        collection = ln.Collection([revised, other], key="lance/release").save()
        with pytest.raises(ValueError, match="separate tables"):
            collection.open()
        assert revised.open_lance().count_rows() == 2
        assert other.open_lance().count_rows() == 1
    finally:
        if collection is not None:
            collection.delete(permanent=True)
        for record in (revised, other, artifact, schema, feature):
            if record is not None:
                record.delete(permanent=True)


def test_append_lance_keeps_one_directory_and_pins_versions():
    artifact = ln.Artifact.from_dataframe(
        _cells(["T cell"]), key="lance/append.lance"
    ).save()
    revised = None
    try:
        revised = artifact.append_lance(_cells(["B cell"]))
        artifact.refresh_from_db()

        assert revised.path == artifact.path
        assert revised.is_latest
        assert not artifact.is_latest
        assert revised.n_observations == 2
        assert artifact.overwrite_versions
        assert revised.overwrite_versions
        assert artifact.open_lance().count_rows() == 1
        assert revised.open_lance().count_rows() == 2
        with pytest.raises(ValueError, match="overwritten"):
            artifact.open()

        again = revised.append_lance(_cells(["NK cell"]))
        revised.refresh_from_db()
        assert again.path == revised.path
        assert again.open_lance().count_rows() == 3
        assert revised.open_lance().count_rows() == 2
        again.delete(permanent=True)
        revised = None
    finally:
        # The latest version owns the shared directory and deletes the family.
        if revised is not None:
            revised.delete(permanent=True)
        elif artifact is not None:
            artifact.delete(permanent=True)
