"""LanceDB table artifacts."""

import shutil
from pathlib import Path
from types import SimpleNamespace
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
    from_path = None
    try:
        assert artifact.suffix == ".lance"
        assert artifact.key == "lance/cells.lance"
        assert artifact.otype == "DataFrame"
        assert artifact.kind == "dataset"
        assert artifact.n_observations == 3
        assert artifact.overwrite_versions is False
        assert artifact.path.is_dir()

        stored = ln.Artifact.get(artifact.uid)
        loaded = cast(pd.DataFrame, stored.load())
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
        shutil.copytree(artifact.path, copied)
        with pytest.raises(ValueError, match="from_lance"):
            ln.Artifact.from_dataframe(copied, key="lance/copied.lance")
        with pytest.raises(ValueError, match=r"\.lance"):
            ln.Artifact.from_lance(df)
        from_path = ln.Artifact.from_lance(
            copied, key="lance/copied.lance", skip_hash_lookup=True
        ).save()
        assert from_path.uid != artifact.uid
        assert from_path.otype == "lance"
        assert from_path.suffix == ".lance"
        assert from_path.n_observations == 3
        table = from_path.open_lance()
        assert table.count_rows() == 3
        from lamindb.core.storage._lance import lance_version

        assert lance_version(from_path) == int(table.version)
        stored_path = ln.Artifact.get(from_path.uid)
        with pytest.raises(ValueError, match="open_lance"):
            stored_path.load()

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
        for record in (from_path, revised, other, artifact, schema, feature):
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
        assert revised.otype == "DataFrame"
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


def test_lance_open_and_append_reject_bad_calls_and_pin_missing_version():
    from lamindb.core.storage._lance import LANCE_VERSION_AUX_KEY, lance_version

    parquet = ln.Artifact.from_dataframe(
        _cells(["parquet"]), key="lance/not-lance.parquet"
    ).save()
    try:
        with pytest.raises(ValueError, match="open_lance"):
            parquet.open_lance()
        with pytest.raises(ValueError, match="append_lance"):
            parquet.append_lance(_cells(["B cell"]))
    finally:
        parquet.delete(permanent=True)

    unsaved = ln.Artifact.from_dataframe(
        _cells(["unsaved-lance"]),
        key="lance/unsaved.lance",
        skip_hash_lookup=True,
    )
    with pytest.raises(ValueError, match="Save the artifact"):
        unsaved.append_lance(_cells(["B cell"]))

    artifact = ln.Artifact.from_dataframe(
        _cells(["unpinned"]), key="lance/unpinned.lance"
    ).save()
    revised = None
    try:
        aux = dict(artifact._aux or {})
        aux.pop(LANCE_VERSION_AUX_KEY, None)
        artifact._aux = aux
        artifact.save()
        assert lance_version(artifact) is None

        revised = artifact.append_lance(_cells(["B cell"]))
        artifact.refresh_from_db()
        pinned = lance_version(artifact)
        assert pinned is not None
        assert pinned != lance_version(revised)
        assert artifact.open_lance().count_rows() == 1
        assert revised.open_lance(version=pinned).count_rows() == 1
        with pytest.raises(ValueError, match="latest"):
            artifact.append_lance(_cells(["NK cell"]))
    finally:
        if revised is not None:
            revised.delete(permanent=True)
        elif artifact is not None:
            artifact.delete(permanent=True)


def test_lance_suffix_from_format_dict():
    from lamindb.core.storage.objects import infer_suffix

    frame = _cells(["T cell"])
    assert infer_suffix(frame, format={"suffix": ".lance"}) == ".lance"
    assert infer_suffix(frame, format=".lance") == ".lance"
    assert infer_suffix(frame, format={"suffix": ".csv"}) == ".csv"


def test_lance_from_polars_frame_and_lazyframe():
    pl = pytest.importorskip("polars")
    frame = pl.DataFrame({"cell_type": ["T cell", None], "n_genes": [1, 2]})
    eager = ln.Artifact.from_dataframe(frame, key="lance/polars.lance").save()
    lazy = ln.Artifact.from_dataframe(frame.lazy(), key="lance/lazy.lance").save()
    try:
        assert eager.n_observations == 2
        assert lazy.n_observations == 2
        assert eager.open_lance().count_rows() == 2
        assert lazy.open_lance().count_rows() == 2
    finally:
        lazy.delete(permanent=True)
        eager.delete(permanent=True)


def test_lance_storage_options_from_s3_filesystem():
    from lamindb.core.storage._lance import lance_storage_options

    anon = ln.UPath(
        "s3://bucket/table.lance?endpoint_url=http://localhost:9000/s3", anon=True
    )
    anon_options = lance_storage_options(anon)
    assert anon_options["aws_endpoint"] == "http://localhost:9000/s3"
    assert anon_options["aws_virtual_hosted_style_request"] == "false"
    assert anon_options["allow_http"] == "true"
    assert anon_options["skip_signature"] == "true"

    keyed = ln.UPath(
        "s3://bucket/table.lance?endpoint_url=https://localhost:9000/s3",
        key="AKIA",
        secret="session-secret",  # noqa: S106
        token="session-token",  # noqa: S106
    )
    keyed_options = lance_storage_options(keyed)
    assert keyed_options["aws_access_key_id"] == "AKIA"
    assert keyed_options["aws_secret_access_key"] == "session-secret"  # noqa: S105
    assert keyed_options["aws_session_token"] == "session-token"  # noqa: S105
    assert keyed_options["aws_endpoint"] == "https://localhost:9000/s3"
    assert "allow_http" not in keyed_options
    assert "skip_signature" not in keyed_options

    assert lance_storage_options(ln.UPath("table.lance")) is None


def test_lance_storage_options_snapshots_refreshable_credentials(monkeypatch):
    from lamindb.core.storage._lance import lance_storage_options

    def provider():
        return {
            "aws_access_key_id": "AKIA",
            "aws_secret_access_key": "secret",
            "aws_session_token": "token",
        }, 123

    monkeypatch.setattr(
        "lamindb.core.storage._polars_lazy_df._polars_options",
        lambda path: {
            "storage_options": {"aws_region": "us-east-1"},
            "credential_provider": provider,
        },
    )
    options = lance_storage_options(SimpleNamespace(protocol="s3"))
    assert options == {
        "aws_region": "us-east-1",
        "aws_access_key_id": "AKIA",
        "aws_secret_access_key": "secret",
        "aws_session_token": "token",
    }


def test_open_lance_table_passes_s3_storage_options(monkeypatch):
    from lamindb.core.storage._lance import open_lance_table

    captured: dict = {}

    class FakeDatabase:
        def open_table(self, name):
            captured["name"] = name
            return "table"

    def connect(uri, **kwargs):
        captured["uri"] = uri
        captured["kwargs"] = kwargs
        return FakeDatabase()

    monkeypatch.setattr(lancedb, "connect", connect)
    monkeypatch.setattr(
        "lamindb.core.storage._lance.lance_storage_options",
        lambda path: {"aws_access_key_id": "AKIA"},
    )
    path = SimpleNamespace(
        parent=SimpleNamespace(as_posix=lambda: "s3://bucket/prefix"),
        stem="cells",
    )
    assert open_lance_table(path) == "table"
    assert captured == {
        "uri": "s3://bucket/prefix",
        "kwargs": {"storage_options": {"aws_access_key_id": "AKIA"}},
        "name": "cells",
    }
