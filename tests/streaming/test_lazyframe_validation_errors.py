from pathlib import Path
from uuid import uuid4

import pytest

pl = pytest.importorskip("polars")


@pytest.fixture
def integer_schema():
    import lamindb as ln

    feature = ln.Feature(
        name=f"stream_policy_{uuid4().hex}", dtype=int, nullable=False
    ).save()
    return ln.Schema(features=[feature], minimal_set=True).save()


def test_keep_preserves_all_rows_after_late_failure(integer_schema, monkeypatch):
    import lamindb as ln
    from lamindb_setup import logger

    name = integer_schema.members.first().name
    frames = [
        pl.LazyFrame({name: [1] * 1000}),
        pl.LazyFrame({name: pl.Series([None], dtype=pl.Int64)}),
        pl.LazyFrame({name: [2] * 1000}),
    ]
    seen_rows = []
    messages = []

    def observe(batch):
        seen_rows.append(batch.height)
        return batch

    first = ln.Artifact.from_lazyframe(
        pl.LazyFrame({name: [0]}),
        key=f"keep-invalid/{uuid4().hex}.parquet",
        schema=integer_schema,
    ).save()
    monkeypatch.setattr(
        logger,
        "warning",
        lambda message, *args, **kwargs: messages.append(str(message)),
    )
    pending = ln.Artifact.from_lazyframe(
        pl.concat(frames, parallel=False).map_batches(
            observe, schema={name: pl.Int64}, streamable=True
        ),
        key=first.key,
        schema=integer_schema,
        on_validation_error="keep",
    )
    assert pending.validation_error is None
    assert not seen_rows
    artifact = pending.save()
    assert sum(seen_rows) == 2001
    assert (
        pl.read_parquet(artifact.path)[name].to_list()
        == [1] * 1000 + [None] + [2] * 1000
    )
    assert artifact.n_observations == 2001
    assert artifact.schema is None
    assert not artifact.schemas.exists()
    assert artifact.stem_uid == first.stem_uid
    first.refresh_from_db()
    assert not first.is_latest
    assert first.schema == integer_schema
    assert pending.validation_error
    artifact.refresh_from_db()
    assert artifact.extra_data["validation"] == {
        "status": "failed",
        "schema_uid": integer_schema.uid,
        "error": pending.validation_error,
    }
    assert sum("Validation failed" in message for message in messages) == 1
    assert pending.save() is artifact
    assert sum(seen_rows) == 2001


def test_keep_preserves_wrong_dtype(integer_schema):
    import lamindb as ln

    name = integer_schema.members.first().name
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({name: ["1", "not-an-int"]}),
        key=f"keep-wrong-dtype/{uuid4().hex}.parquet",
        schema=integer_schema,
        on_validation_error="keep",
    )
    assert not pending.path.exists()
    artifact = pending.save()
    assert pending.validation_error
    assert pl.read_parquet(artifact.path).to_dict(as_series=False) == {
        name: ["1", "not-an-int"]
    }
    assert artifact.schema is None
    assert artifact.extra_data["validation"]["status"] == "failed"


def test_keep_handles_missing_columns_during_planning(integer_schema):
    import lamindb as ln

    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({"unexpected": [1]}),
        key=f"keep-missing-column/{uuid4().hex}.parquet",
        schema=integer_schema,
        on_validation_error="keep",
    )
    assert pending.validation_error
    assert not pending.path.exists()
    artifact = pending.save()
    assert pl.read_parquet(artifact.path).to_dict(as_series=False) == {
        "unexpected": [1]
    }
    assert artifact.schema is None
    assert artifact.extra_data["validation"]["status"] == "failed"


def test_keep_categorical_failure_does_not_attach_partial_annotations():
    import lamindb as ln

    label = ln.ULabel(name=f"policy_label_{uuid4().hex}").save()
    name = f"policy_category_{uuid4().hex}"
    schema = ln.Schema(
        features=[ln.Feature(name=name, dtype="cat[ULabel]").save()]
    ).save()
    invalid = f"unknown_label_{uuid4().hex}"
    pending = ln.Artifact.from_lazyframe(
        pl.concat(
            [
                pl.LazyFrame({name: [label.name] * 1000}),
                pl.LazyFrame({name: [invalid]}),
            ],
            parallel=False,
        ),
        key=f"keep-invalid-category/{uuid4().hex}.parquet",
        schema=schema,
        on_validation_error="keep",
    )
    artifact = pending.save()
    assert pending.validation_error
    assert artifact.schema is None
    assert not artifact.ulabels.exists()
    assert not artifact.schemas.exists()
    assert pl.read_parquet(artifact.path)[name].to_list() == [label.name] * 1000 + [
        invalid
    ]


@pytest.mark.parametrize("coercion", ["schema", "feature"])
def test_keep_rejects_coercion_before_execution_or_reservation(coercion):
    import lamindb as ln
    from lamindb.models._artifact_uid_reservation import ArtifactUIDReservation

    name = f"coercing_policy_{uuid4().hex}"
    schema = ln.Schema(
        features=[
            ln.Feature(name=name, dtype=int, coerce=coercion == "feature").save()
        ],
        coerce=coercion == "schema",
    ).save()
    seen_rows = []

    def observe(batch):
        seen_rows.append(batch.height)
        return batch

    frame = pl.LazyFrame({name: ["1", "invalid"]}).map_batches(
        observe, schema={name: pl.String}, streamable=True
    )
    count = ArtifactUIDReservation.objects.count()
    with pytest.raises(ValueError, match="non-coercing schema"):
        ln.Artifact.from_lazyframe(
            frame,
            key=f"reject-coercion/{uuid4().hex}.parquet",
            schema=schema,
            on_validation_error="keep",
        )
    assert not seen_rows
    assert ArtifactUIDReservation.objects.count() == count


def test_raise_cleans_partial_output_and_preserves_head(integer_schema, monkeypatch):
    import lamindb as ln
    from lamindb.errors import ValidationError
    from lamindb.models._artifact_uid_reservation import ArtifactUIDReservation

    name = integer_schema.members.first().name
    first = ln.Artifact.from_lazyframe(
        pl.LazyFrame({name: [1]}),
        key=f"validation-head/{uuid4().hex}.parquet",
        schema=integer_schema,
    ).save()
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({name: pl.Series([None], dtype=pl.Int64)}),
        key=first.key,
        schema=integer_schema,
        on_validation_error="raise",
    )

    def sink_partial(frame, path, **kwargs):
        Path(path).write_bytes(b"partial parquet")
        frame.collect(engine="streaming")
        pytest.fail("The invalid batch must abort the sink")

    monkeypatch.setattr(pl.LazyFrame, "sink_parquet", sink_partial)
    with pytest.raises(ValidationError):
        pending.save()
    assert not pending.path.exists()
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    assert (
        ArtifactUIDReservation.objects.get(uid=pending.path.stem).state == "discarded"
    )
    first.refresh_from_db()
    assert first.is_latest
    assert first.schema == integer_schema


def test_cleanup_failure_preserves_validation_error_and_can_retry(
    integer_schema, monkeypatch
):
    import lamindb as ln
    from lamindb.errors import ValidationError
    from lamindb.models._artifact_uid_reservation import ArtifactUIDReservation

    name = integer_schema.members.first().name
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({name: pl.Series([None], dtype=pl.Int64)}),
        key=f"validation-cleanup/{uuid4().hex}.parquet",
        schema=integer_schema,
    )

    def fail_cleanup(*args, **kwargs):
        raise OSError("cleanup unavailable")

    with monkeypatch.context() as patch:
        patch.setattr(type(pending.path), "unlink", fail_cleanup)
        with pytest.raises(ValidationError) as error:
            pending.save()
    assert isinstance(error.value.__cause__, OSError)
    assert "cleanup unavailable" in str(error.value.__cause__)
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    assert ArtifactUIDReservation.objects.get(uid=pending.path.stem).state == "failed"
    pending.discard()
    assert not pending.path.exists()
    assert (
        ArtifactUIDReservation.objects.get(uid=pending.path.stem).state == "discarded"
    )


@pytest.mark.parametrize("on_validation_error", ["raise", "keep"])
def test_failure_policy_does_not_swallow_io_errors(
    integer_schema, monkeypatch, on_validation_error
):
    import lamindb as ln

    name = integer_schema.members.first().name
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({name: [1]}),
        key=f"validation-io/{uuid4().hex}.parquet",
        schema=integer_schema,
        on_validation_error=on_validation_error,
    )

    def fail_sink(frame, path, **kwargs):
        Path(path).write_bytes(b"partial parquet")
        raise OSError("upload unavailable")

    monkeypatch.setattr(pl.LazyFrame, "sink_parquet", fail_sink)
    with pytest.raises(OSError, match="upload unavailable"):
        pending.save()
    assert pending.path.exists()
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()
    assert pending.validation_error is None
    pending.discard()


def test_keep_deduplication_does_not_modify_existing_validation(integer_schema):
    import lamindb as ln

    name = integer_schema.members.first().name
    frame = pl.LazyFrame({name: ["1", "invalid"]})
    first = ln.Artifact.from_lazyframe(
        frame, key=f"keep-duplicate/{uuid4().hex}.parquet"
    ).save()
    first.extra_data = {"original": True}
    first.save()
    pending = ln.Artifact.from_lazyframe(
        frame,
        key=first.key,
        schema=integer_schema,
        on_validation_error="keep",
    )
    assert pending.save().uid == first.uid
    assert not pending.path.exists()
    assert pending.validation_error
    first.refresh_from_db()
    assert first.extra_data == {"original": True}
    assert first.schema is None
    assert first.is_latest


def test_keep_validation_failure_survives_database_retry(integer_schema, monkeypatch):
    import lamindb as ln
    from lamindb.models._artifact_uid_reservation import ArtifactUIDReservation

    name = integer_schema.members.first().name
    pending = ln.Artifact.from_lazyframe(
        pl.LazyFrame({name: pl.Series([1, None], dtype=pl.Int64)}),
        key=f"keep-retry/{uuid4().hex}.parquet",
        schema=integer_schema,
        on_validation_error="keep",
    )
    save_reservation = ArtifactUIDReservation.save

    def fail_finalization(reservation, *args, **kwargs):
        if reservation.state == "complete":
            raise OSError("database unavailable")
        return save_reservation(reservation, *args, **kwargs)

    with monkeypatch.context() as patch:
        patch.setattr(ArtifactUIDReservation, "save", fail_finalization)
        with pytest.raises(OSError, match="database unavailable"):
            pending.save()
    assert not ln.Artifact.objects.filter(uid=pending.path.stem).exists()

    def unexpected_sink(*args, **kwargs):
        pytest.fail("A completed upload must not be repeated")

    monkeypatch.setattr(pl.LazyFrame, "sink_parquet", unexpected_sink)
    artifact = pending.save()
    assert artifact.schema is None
    assert artifact.extra_data["validation"]["error"] == pending.validation_error
