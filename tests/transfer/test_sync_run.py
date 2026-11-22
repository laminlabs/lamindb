from __future__ import annotations

import tempfile
import time
from pathlib import Path

import lamindb as ln
import pytest


def _source() -> str:
    return f"{ln.setup.settings.user.handle}/testdb1"


def _text_file(name: str, content: str) -> Path:
    path = Path(tempfile.gettempdir()) / name
    path.write_text(content)
    return path


def _sync_runs(uid: str):
    return ln.Run.filter(entrypoint="sync", params__uid=uid)


def test_sync_creates_a_run_per_call_and_save_reuses_one():
    unique = time.time_ns()
    ln.connect("testdb1")
    first = ln.Record(name=f"sync-run-a-{unique}", is_type=True).save()
    second = ln.Record(name=f"sync-run-b-{unique}", is_type=True).save()
    third = ln.Record(name=f"sync-run-c-{unique}", is_type=True).save()
    fourth = ln.Record(name=f"sync-run-d-{unique}", is_type=True).save()
    artifact = ln.Artifact(
        _text_file(f"sync-map-{unique}.txt", f"sync map {unique}\n"),
        key=f"sync-map-{unique}.txt",
    ).save()
    ln.connect("testdb2")
    source = _source()

    synced_first = ln.core.sync(registry=ln.Record, uid=first.uid, source_db=source)
    synced_second = ln.core.sync(registry=ln.Record, uid=second.uid, source_db=source)

    assert synced_first.run.uid != synced_second.run.uid
    assert synced_first.run.transform.uid == synced_second.run.transform.uid
    assert synced_first.run.transform.key == (
        f"__lamindb_transfer__/{synced_first.run.transform.uid[:-4]}"
    )
    assert synced_first.run.entrypoint == "sync"
    assert synced_first.run.initiated_by_run is None
    assert synced_first.run.status == "completed"
    assert synced_first.run.finished_at is not None
    assert synced_first.run.params == {
        "registry": "Record",
        "uid": first.uid,
        "source_db": source,
        "transfer": None,
        "depth": 0,
    }
    assert ln.context.run is None

    again = ln.core.sync(registry=ln.Record, uid=first.uid, source_db=source)
    assert again.id == synced_first.id
    # A Record re-sync fills the row, so `.run` is this call's run.
    assert again.run.uid != synced_first.run.uid
    assert again.run.entrypoint == "sync"
    assert _sync_runs(first.uid).count() == 2

    synced_artifact = ln.core.sync(
        registry=ln.Artifact, uid=artifact.uid, source_db=source
    )
    again_artifact = ln.core.sync(
        registry=ln.Artifact, uid=artifact.uid, source_db=source
    )
    assert again_artifact.id == synced_artifact.id
    assert again_artifact.run.uid == synced_artifact.run.uid
    assert _sync_runs(artifact.uid).count() == 2

    db = ln.DB(source)
    saved_third = db.Record.get(third.uid).save()
    saved_fourth = db.Record.get(fourth.uid).save()
    assert saved_third.run.uid == saved_fourth.run.uid
    assert saved_third.run.entrypoint is None
    assert saved_third.run.uid != synced_first.run.uid
    assert saved_third.run.initiated_by_run is None
    assert saved_third.run.status == "completed"


def test_sync_depth_and_collection_share_one_run():
    unique = time.time_ns()
    ln.connect("testdb1")
    parent = ln.Record(name=f"sync-depth-parent-{unique}", is_type=True).save()
    child = ln.Record(
        name=f"sync-depth-child-{unique}", is_type=True, type=parent
    ).save()
    row = ln.Record(name=f"sync-depth-row-{unique}", type=child).save()
    artifact = ln.Artifact(
        _text_file(f"sync-collection-{unique}.txt", f"sync collection {unique}\n"),
        key=f"sync-collection-{unique}.txt",
    ).save()
    collection = ln.Collection(artifact, key=f"sync-collection-{unique}").save()
    ln.connect("testdb2")
    source = _source()

    synced = ln.core.sync(registry=ln.Record, uid=parent.uid, source_db=source, depth=2)
    assert _sync_runs(parent.uid).count() == 1
    assert ln.Record.get(child.uid).run.uid == synced.run.uid
    assert ln.Record.get(row.uid).run.uid == synced.run.uid
    assert synced.run.params["depth"] == 2

    synced_collection = ln.core.sync(
        registry=ln.Collection, uid=collection.uid, source_db=source
    )
    member = ln.Artifact.get(uid=artifact.uid)
    assert member.run.uid == synced_collection.run.uid
    assert _sync_runs(collection.uid).count() == 1
    assert synced_collection.run.uid != synced.run.uid


def test_sync_run_is_initiated_by_the_tracked_run():
    unique = time.time_ns()
    ln.connect("testdb1")
    record = ln.Record(name=f"sync-tracked-{unique}", is_type=True).save()
    other = ln.Record(name=f"sync-tracked-save-{unique}", is_type=True).save()
    ln.connect("testdb2")
    transform = ln.Transform(key=f"sync-tracked-{unique}.py").save()
    source = _source()
    try:
        ln.track(transform=transform)
        parent = ln.context.run
        synced = ln.core.sync(registry=ln.Record, uid=record.uid, source_db=source)
        assert ln.context.run.uid == parent.uid
        assert synced.run.initiated_by_run.uid == parent.uid
        assert synced.run.entrypoint == "sync"
        assert synced.run.transform.uid != parent.transform.uid

        saved = ln.DB(source).Record.get(other.uid).save()
        assert saved.run.initiated_by_run.uid == parent.uid
        assert saved.run.entrypoint is None
        assert saved.run.uid != synced.run.uid
        ln.finish()
    finally:
        ln.context._run = None


def test_sync_failure_marks_the_run_errored(monkeypatch: pytest.MonkeyPatch):
    unique = time.time_ns()
    ln.connect("testdb1")
    record = ln.Record(name=f"sync-fails-{unique}", is_type=True).save()
    ln.connect("testdb2")
    original = ln.Record.save

    def fail_transfer(self, *args, **kwargs):
        db = getattr(getattr(self, "_state", None), "db", None)
        if db not in (None, "default"):
            raise RuntimeError("sync failed")
        return original(self, *args, **kwargs)

    monkeypatch.setattr(ln.Record, "save", fail_transfer)
    with pytest.raises(RuntimeError, match="sync failed"):
        ln.core.sync(registry=ln.Record, uid=record.uid, source_db=_source())

    assert ln.Record.filter(uid=record.uid).one_or_none() is None
    run = _sync_runs(record.uid).one()
    assert run.status == "errored"
    assert run.finished_at is not None
    assert run.entrypoint == "sync"
