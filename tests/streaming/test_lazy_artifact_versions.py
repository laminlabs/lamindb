from uuid import uuid4

import pytest

pl = pytest.importorskip("polars")


def lazy(key=None, **kwargs):
    import lamindb as ln

    return ln.Artifact.from_lazy(
        suffix=".parquet",
        overwrite_versions=False,
        key=key or f"streaming/{uuid4().hex}.parquet",
        deduplicate="key",
        **kwargs,
    )


def write(pending, values=(1, 2)):
    pl.LazyFrame({"value": values}).sink_parquet(pending.path.as_posix())
    return pending


def test_managed_path_and_no_cache(monkeypatch):
    import lamindb as ln
    from django.db import connection
    from lamindb.models.artifact import _s

    def unexpected_cache(*args, **kwargs):
        pytest.fail("Streaming must not serialize to the local cache")

    monkeypatch.setattr(_s(), "write_to_disk", unexpected_cache)
    pending = lazy()
    managed_path = pending.path
    assert not managed_path.exists()
    assert not connection.in_atomic_block
    artifact = write(pending).save()
    assert artifact.path == managed_path
    assert artifact.path.name == f"{artifact.uid}.parquet"
    assert artifact._real_key is None
    assert artifact._key_is_virtual
    assert not artifact._storage_ongoing
    assert artifact.n_observations == 2
    assert artifact.size == artifact.path.stat().st_size
    assert artifact.hash is not None
    assert pl.read_parquet(artifact.path).to_dict(as_series=False) == {"value": [1, 2]}
    assert ln.Artifact.get(path=artifact.path) == artifact
    assert pending.save() is artifact
    artifact.key = f"renamed/{uuid4().hex}.parquet"
    artifact.save()
    assert artifact.path == managed_path
    with pytest.raises(ValueError, match="published"):
        pending.discard()


def test_same_key_dedup_and_different_key_duplicates():
    first = write(lazy()).save()
    duplicate = write(lazy(key=first.key))
    assert duplicate.save().uid == first.uid
    assert duplicate.path != first.path
    assert not duplicate.path.exists()
    assert write(lazy()).save().uid != first.uid


def test_abandoned_writes_do_not_reserve_uids_or_insert_rows():
    import lamindb as ln

    first = write(lazy()).save()
    count = ln.Artifact.objects.count()
    abandoned = lazy(key=first.key)
    abandoned_uid = abandoned.path.stem
    assert ln.Artifact.objects.count() == count
    abandoned.discard()
    second = write(lazy(key=first.key), [3]).save()
    assert second.stem_uid == first.stem_uid
    assert second.uid == abandoned_uid
    assert second.uid[-4:] == "0001"
    first.refresh_from_db()
    assert not first.is_latest
    assert second.is_latest


def test_normal_creation_is_not_changed_by_pending_stream():
    import lamindb as ln

    first = write(lazy()).save()
    pending = lazy(key=first.key)
    normal = ln.Artifact.from_dataframe(
        pl.LazyFrame({"value": [7]}), key=first.key, skip_hash_lookup=True
    ).save()
    assert normal.stem_uid == first.stem_uid
    assert normal.uid == pending.path.stem
    assert normal.uid[-4:] == "0001"
    with pytest.raises(ValueError, match="published"):
        pending.discard()


def test_pending_stream_cannot_overwrite_or_discard_published_uid():
    import lamindb as ln
    from lamindb.errors import IntegrityError

    first = write(lazy()).save()
    normal = ln.Artifact.from_dataframe(
        pl.LazyFrame({"value": [7]}), key=first.key, skip_hash_lookup=True
    )
    pending = lazy(key=first.key)
    assert normal.uid == pending.path.stem
    normal.save()
    with pytest.raises(IntegrityError, match="published by another writer"):
        pending.save()
    with pytest.raises(ValueError, match="published"):
        pending.discard()
    assert pl.read_parquet(normal.path)["value"].to_list() == [7]


def test_lazyframe_refuses_occupied_path_before_sinking():
    import lamindb as ln
    from lamindb.errors import IntegrityError

    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [1]}), key=f"occupied/{uuid4().hex}.parquet"
    )
    pending.path.write_bytes(b"unfinished data")
    with pytest.raises(IntegrityError, match="refusing to overwrite"):
        pending.save()
    assert pending.path.read_bytes() == b"unfinished data"
    pending.discard()


def test_existing_partial_revision_must_be_cleaned_up():
    from lamindb.errors import IntegrityError

    first = write(lazy()).save()
    pending = write(lazy(key=first.key), [8])
    with pytest.raises(IntegrityError, match="destination already exists"):
        lazy(key=first.key)
    pending.discard()


def test_two_pending_revisions_have_no_reservation():
    first = write(lazy()).save()
    left = lazy(key=first.key)
    right = lazy(key=first.key)
    assert left.path == right.path
    left.discard()
    right.discard()


def test_explicit_revises_works_when_writes_are_serialized():
    import lamindb as ln

    first = write(lazy()).save()
    pending = ln.Artifact.from_lazyframe(pl.LazyFrame({"value": [7]}), revises=first)
    assert pending.path.stem == f"{first.stem_uid}0001"
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    revised = pending.save()
    assert revised.uid == pending.path.stem
    assert revised.key == first.key
    assert pl.read_parquet(first.path)["value"].to_list() == [1, 2]
    assert pl.read_parquet(revised.path)["value"].to_list() == [7]
    first.refresh_from_db()
    assert not first.is_latest
    assert revised.is_latest


def test_duplicate_published_during_streaming():
    pending = lazy()
    concurrent = write(lazy(key=pending.kwargs["key"]), [9]).save()
    result = write(pending, [9]).save()
    assert result.uid == concurrent.uid
    assert result.uid != pending.path.stem
    concurrent.refresh_from_db()
    assert concurrent.is_latest
    assert not pending.path.exists()


def test_different_version_published_during_streaming():
    import lamindb as ln
    from lamindb.errors import IntegrityError

    pending = write(lazy(), [9])
    concurrent = write(lazy(key=pending.kwargs["key"]), [8]).save()
    with pytest.raises(IntegrityError, match="different version"):
        pending.save()
    assert pending.path.exists()
    concurrent.refresh_from_db()
    assert concurrent.is_latest
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    pending.discard()
    assert not pending.path.exists()


def test_duplicate_of_previous_version_does_not_change_head():
    first = write(lazy()).save()
    current = write(lazy(key=first.key), [8]).save()
    pending = write(lazy(key=first.key))
    assert pending.save().uid == first.uid
    first.refresh_from_db()
    current.refresh_from_db()
    assert not first.is_latest
    assert current.is_latest


def test_invalid_or_incomplete_upload_is_not_published():
    import lamindb as ln

    pending = lazy()
    pending.path.write_bytes(b"partial")
    with pytest.raises(Exception, match="Parquet"):
        pending.save()
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    write(pending)
    assert pending.save().uid == pending.path.stem


def test_discard_removes_partial_upload():
    pending = lazy()
    pending.path.write_bytes(b"partial")
    pending.discard()
    assert not pending.path.exists()
    with pytest.raises(ValueError, match="discarded"):
        pending.save()


def test_metadata_failure_can_retry_without_rewriting(monkeypatch):
    import pyarrow.parquet as pq

    pending = write(lazy())
    read_metadata = pq.read_metadata
    with monkeypatch.context() as patch:

        def fail(*args, **kwargs):
            raise OSError("footer unavailable")

        patch.setattr(pq, "read_metadata", fail)
        with pytest.raises(OSError, match="footer unavailable"):
            pending.save()
    assert pq.read_metadata is read_metadata
    assert pending.save().path == pending.path


def test_cleanup_failure_can_retry(monkeypatch):
    from upath import UPath

    first = write(lazy()).save()
    pending = write(lazy(key=first.key))
    with monkeypatch.context() as patch:

        def fail(*args, **kwargs):
            raise OSError("cleanup unavailable")

        patch.setattr(type(UPath(pending.path)), "unlink", fail)
        with pytest.raises(OSError, match="cleanup unavailable"):
            pending.save()
    assert pending.path.exists()
    assert pending.save().uid == first.uid
    assert not pending.path.exists()


def test_legacy_from_lazy_unchanged():
    import lamindb as ln

    pending = ln.Artifact.from_lazy(
        suffix=".parquet",
        overwrite_versions=False,
        key=f"legacy/{uuid4().hex}.parquet",
    )
    pending.path.parent.mkdir(parents=True, exist_ok=True)
    artifact = write(pending).save()
    assert artifact.uid == pending.path.stem
    with pytest.raises(ValueError, match="deduplicate"):
        pending.discard()


def test_missing_file_can_be_written_then_saved():
    pending = lazy()
    with pytest.raises(FileNotFoundError):
        pending.save()
    assert write(pending).save().path == pending.path


def test_streaming_rejects_outer_transaction():
    from django.db import transaction

    with transaction.atomic():
        with pytest.raises(ValueError, match="outside"):
            lazy()
    pending = write(lazy())
    with transaction.atomic():
        with pytest.raises(ValueError, match="outside"):
            pending.save()
    pending.discard()


def test_validates_and_coerces_inside_single_source_execution():
    import lamindb as ln

    name = f"stream_int_{uuid4().hex}"
    schema = ln.Schema(
        features=[ln.Feature(name=name, dtype=int, coerce=True).save()],
        maximal_set=True,
    ).save()
    seen_rows = []

    def observe(batch):
        seen_rows.append(batch.height)
        return batch

    frame = pl.LazyFrame({name: ["1", "2"] * 50000}).map_batches(
        observe, schema={name: pl.String}, streamable=True
    )
    pending = ln.Artifact.from_lazyframe(
        frame,
        key=f"validated/{uuid4().hex}.parquet",
        schema=schema,
    )
    assert not seen_rows
    artifact = pending.save()
    assert sum(seen_rows) == 100000
    saved = pl.read_parquet(artifact.path)
    assert saved.schema[name] == pl.Int64
    assert saved[name].to_list() == [1, 2] * 50000
    assert artifact.schema == schema
    assert artifact.n_observations == 100000
    pending.save()
    assert sum(seen_rows) == 100000


def test_invalid_batch_prevents_registration():
    import lamindb as ln
    from lamindb.errors import ValidationError

    name = f"invalid_stream_{uuid4().hex}"
    schema = ln.Schema(
        features=[ln.Feature(name=name, dtype=int, coerce=True).save()]
    ).save()
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({name: ["1", "not-an-int"]}),
        key=f"invalid/{uuid4().hex}.parquet",
        schema=schema,
    )
    with pytest.raises(ValidationError):
        pending.save()
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    assert not pending.path.exists()
    assert pending.validation_error is not None
    with pytest.raises(ValueError, match="discarded"):
        pending.save()
    pending.discard()


def test_empty_frame_preserves_validated_schema():
    import lamindb as ln

    name = f"empty_stream_{uuid4().hex}"
    schema = ln.Schema(
        features=[ln.Feature(name=name, dtype=int, coerce=True).save()]
    ).save()
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame(schema={name: pl.String}),
        key=f"empty/{uuid4().hex}.parquet",
        schema=schema,
    )
    artifact = pending.save()
    assert artifact.schema == schema
    assert artifact.n_observations == 0
    assert pl.read_parquet(artifact.path).schema[name] == pl.Int64


@pytest.mark.parametrize("on_validation_error", ["raise", "keep"])
def test_frame_without_schema_sinks_once(on_validation_error):
    import lamindb as ln

    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [1, 2]}),
        key=f"unvalidated/{uuid4().hex}.parquet",
        on_validation_error=on_validation_error,
    )
    artifact = pending.save()
    assert artifact.path == pending.path
    assert artifact.n_observations == 2
    assert artifact.schema is None


def test_lazyframe_defaults_and_parquet_options(monkeypatch):
    import lamindb as ln
    import pyarrow.parquet as pq
    from lamindb.models.artifact import _s

    def unexpected_cache(*args, **kwargs):
        pytest.fail("Streaming must not serialize to the local cache")

    monkeypatch.setattr(_s(), "write_to_disk", unexpected_cache)
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [1, 2]}),
        key=f"lazyframe-options/{uuid4().hex}.parquet",
        parquet_kwargs={"compression": "uncompressed"},
    )
    assert not pending.path.exists()
    artifact = pending.save()
    assert artifact.suffix == ".parquet"
    assert artifact.path.name == f"{artifact.uid}.parquet"
    assert (
        pq.read_metadata(artifact.path).row_group(0).column(0).compression
        == "UNCOMPRESSED"
    )


def test_keyless_lazyframes_have_independent_families():
    import lamindb as ln

    first = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [101]}), description="First keyless stream"
    ).save()
    second = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [102]}), description="Second keyless stream"
    ).save()
    duplicate = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [101]}), description="Duplicate keyless stream"
    )
    assert duplicate.save() == first
    assert not duplicate.path.exists()
    first.refresh_from_db()
    assert first.is_latest and second.is_latest
    assert first.stem_uid != second.stem_uid
    assert first.key is None and second.key is None

    revised = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [103]}), revises=first
    ).save()
    assert revised.stem_uid == first.stem_uid
    first.refresh_from_db()
    second.refresh_from_db()
    assert not first.is_latest
    assert second.is_latest and revised.is_latest


def test_lazyframe_inherits_revision_key_and_external_features():
    import lamindb as ln

    first = write(lazy()).save()
    feature = ln.Feature(name=f"stream_external_{uuid4().hex}", dtype=int).save()
    revised = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [42]}),
        revises=first,
        features={feature.name: 7},
    ).save()
    assert revised.key == first.key
    assert revised.stem_uid == first.stem_uid
    assert revised.features.get_values()[feature.name] == 7


def test_lazyframe_valid_features_matches_dataframe_restrictions():
    import lamindb as ln
    from lamindb.errors import InvalidArgument

    name = f"stream_valid_feature_{uuid4().hex}"
    ln.Feature(name=name, dtype=int).save()
    for constructor in (ln.Artifact.from_dataframe, ln.Artifact.from_lazyframe):
        with pytest.warns(DeprecationWarning):
            with pytest.raises(InvalidArgument, match="flexible schema"):
                constructor(
                    pl.LazyFrame({name: [1, 2]}),
                    key=f"valid-features/{uuid4().hex}.parquet",
                    schema="valid_features",
                )


@pytest.mark.parametrize("data", [pl.DataFrame({"value": [1]}), "input.parquet", None])
def test_lazyframe_rejects_non_lazyframes(data):
    import lamindb as ln

    with pytest.raises(ValueError, match="requires a Polars LazyFrame"):
        ln.Artifact.from_lazyframe(data, key="invalid.parquet")


@pytest.mark.parametrize(
    "options",
    [
        {"key": "file.csv"},
        {"key": "file.parquet", "format": ".csv"},
        {"key": "file.parquet", "format": {"suffix": ".csv"}},
        {"key": "file.parquet", "format": ".zarr"},
        {"key": "file.arrow"},
        {"key": "file.parquet", "overwrite_versions": True},
        {"key": "file.parquet", "key_is_virtual": False},
        {"key": "file.parquet", "csv_kwargs": {"separator": ";"}},
        {"key": "file.parquet", "parquet_kwargs": {"lazy": True}},
        {"key": "file.parquet", "parquet_kwargs": {"path": "other.parquet"}},
        {"key": "file.parquet", "on_validation_error": "ignore"},
    ],
)
def test_lazyframe_rejects_unsupported_output_before_inserting(options):
    import lamindb as ln

    count = ln.Artifact.objects.count()
    with pytest.raises(ValueError):
        ln.Artifact.from_lazyframe(pl.LazyFrame({"value": [1]}), **options)
    assert ln.Artifact.objects.count() == count


def test_from_lazy_does_not_accept_a_frame():
    import lamindb as ln

    frame = pl.LazyFrame({"value": [1]})
    with pytest.raises(TypeError, match="use from_lazyframe"):
        ln.Artifact.from_lazy(frame, False)
    with pytest.raises(TypeError, match="use from_lazyframe"):
        ln.Artifact.from_lazy(".parquet", False, frame=frame)
    with pytest.raises(TypeError):
        ln.Artifact.from_lazyframe(frame, key="file.parquet", suffix=".parquet")


@pytest.mark.parametrize("on_validation_error", ["raise", "keep"])
def test_categorical_annotations_include_all_batches(on_validation_error):
    import lamindb as ln

    labels = [
        ln.ULabel(name=f"stream_label_{uuid4().hex}").save(),
        ln.ULabel(name=f"stream_label_{uuid4().hex}").save(),
    ]
    name = f"stream_category_{uuid4().hex}"
    schema = ln.Schema(
        features=[ln.Feature(name=name, dtype="cat[ULabel]").save()]
    ).save()
    seen_batches = []

    def observe(batch):
        seen_batches.append(batch.height)
        return batch

    frame = pl.concat(
        [pl.LazyFrame({name: [label.name] * 1000}) for label in labels],
        parallel=False,
    ).map_batches(observe, schema={name: pl.String}, streamable=True)
    artifact = ln.Artifact.from_lazyframe(
        frame,
        key=f"categorical/{uuid4().hex}.parquet",
        schema=schema,
        on_validation_error=on_validation_error,
    ).save()
    assert len(seen_batches) >= 2
    assert sum(seen_batches) == 2000
    assert set(artifact.ulabels.values_list("uid", flat=True)) == {
        label.uid for label in labels
    }
    assert artifact.schema == schema
    assert artifact.extra_data is None


def test_invalid_late_batch_is_not_published():
    import lamindb as ln
    from lamindb.errors import ValidationError

    name = f"late_invalid_{uuid4().hex}"
    schema = ln.Schema(
        features=[ln.Feature(name=name, dtype=int, coerce=True).save()]
    ).save()
    frame = pl.concat(
        [pl.LazyFrame({name: ["1"] * 1000}), pl.LazyFrame({name: ["invalid"]})],
        parallel=False,
    )
    pending = ln.Artifact.from_lazyframe(
        frame,
        key=f"late-invalid/{uuid4().hex}.parquet",
        schema=schema,
    )
    with pytest.raises(ValidationError):
        pending.save()
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    assert not pending.path.exists()
    pending.discard()


def test_database_failure_retries_without_resinking(monkeypatch):
    import lamindb as ln

    first = write(lazy()).save()
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"value": [7]}),
        key=first.key,
    )
    with monkeypatch.context() as patch:
        save_artifact = ln.Artifact.save

        def fail_finalization(artifact, *args, **kwargs):
            result = save_artifact(artifact, *args, **kwargs)
            if artifact.uid == pending.path.stem:
                raise OSError("database unavailable")
            return result

        patch.setattr(ln.Artifact, "save", fail_finalization)
        with pytest.raises(OSError, match="database unavailable"):
            pending.save()
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    first.refresh_from_db()
    assert first.is_latest

    def unexpected_sink(*args, **kwargs):
        pytest.fail("A completed upload must not be repeated")

    monkeypatch.setattr(pl.LazyFrame, "sink_parquet", unexpected_sink)
    result = pending.save()
    assert result.uid == pending.path.stem
    first.refresh_from_db()
    assert not first.is_latest
    assert result.is_latest


def test_branch_heads_remain_independent():
    import lamindb as ln

    main = write(lazy()).save()
    branch = ln.Branch(name=f"streaming-branch-{uuid4().hex}").save()
    alternate = write(lazy(key=main.key, branch=branch), [7]).save()
    updated = write(lazy(key=main.key), [8]).save()
    main.refresh_from_db()
    alternate.refresh_from_db()
    assert not main.is_latest
    assert alternate.is_latest
    assert updated.is_latest
    assert main.stem_uid == alternate.stem_uid == updated.stem_uid
    assert updated.uid[-4:] == "0002"


def test_legacy_folder_flow_unchanged():
    import lamindb as ln

    pending = ln.Artifact.from_lazy(
        suffix="",
        overwrite_versions=True,
        key=f"legacy-folder-{uuid4().hex}",
    )
    pending.path.mkdir(parents=True)
    (pending.path / "data.txt").write_text("folder content")
    artifact = pending.save()
    assert artifact.n_files == 1
    assert artifact.uid[:16] == pending.path.name
    assert artifact.path == pending.path


def test_legacy_positional_arguments():
    import lamindb as ln

    pending = ln.Artifact.from_lazy(
        ".parquet", False, f"legacy-positional/{uuid4().hex}.parquet"
    )
    pending.path.parent.mkdir(parents=True, exist_ok=True)
    assert write(pending).save().path == pending.path


def test_schema_without_frame_rejected():
    import lamindb as ln

    schema = ln.Schema(
        features=[ln.Feature(name=f"unused_{uuid4().hex}", dtype=int).save()]
    ).save()
    with pytest.raises(ValueError, match="from_lazyframe"):
        ln.Artifact.from_lazy(
            suffix=".parquet",
            overwrite_versions=False,
            key="no-frame.parquet",
            schema=schema,
        )


@pytest.mark.parametrize(
    "options",
    [
        {"suffix": ".zarr", "overwrite_versions": True, "key": "folder.zarr"},
        {"suffix": ".parquet", "overwrite_versions": True, "key": "file.parquet"},
        {"suffix": ".csv", "overwrite_versions": False, "key": "file.csv"},
        {"suffix": ".parquet", "overwrite_versions": False},
        {
            "suffix": ".parquet",
            "overwrite_versions": False,
            "key": "file.parquet",
            "key_is_virtual": False,
        },
    ],
)
def test_unsupported_output_fails_before_upload(options):
    import lamindb as ln

    with pytest.raises(ValueError):
        ln.Artifact.from_lazy(deduplicate="key", **options)
