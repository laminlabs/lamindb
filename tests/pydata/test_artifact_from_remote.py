import json
from pathlib import Path
from unittest.mock import patch
from uuid import uuid4

import lamindb as ln
import pandas as pd
import pytest
from lamindb_setup.core.upath import UPath


@pytest.fixture
def target_storage(tmp_path):
    root = tmp_path / "from-remote-storage"
    root.mkdir()
    storage = ln.Storage(root=root.as_posix(), type="local").save()
    yield storage
    for artifact in ln.Artifact.objects.filter(storage=storage).order_by("-id"):
        artifact.delete(permanent=True, storage=True)
    storage.delete()


@pytest.fixture
def source():
    root = UPath(f"memory://from-remote-{uuid4().hex}")
    (root / "file.txt").write_bytes(uuid4().hex.encode())
    (root / "folder" / "a.txt").write_bytes(uuid4().hex.encode())
    (root / "folder" / "sub" / "b.txt").write_bytes(b"bb")
    yield root
    root.fs.rm(root.as_posix(), recursive=True)


@pytest.mark.parametrize(
    "is_folder,bulk", [(False, False), (False, True), (True, False)]
)
def test_from_remote_copies_into_storage(source, target_storage, is_folder, bulk):
    path = source / ("folder" if is_folder else "file.txt")
    artifact = ln.Artifact.from_remote(
        path, key=f"imports/{path.name}", storage=target_storage
    )
    assert artifact.path == path
    source_stats = artifact.extra_data["source"]
    assert source_stats["url"] == path.as_posix()
    assert source_stats["size"] == artifact.size
    if bulk:
        ln.save([artifact])
    else:
        artifact.save()

    loaded = ln.Artifact.get(artifact.uid)
    assert loaded.path.is_relative_to(UPath(target_storage.root))
    assert loaded.extra_data["source"] == source_stats
    assert loaded._hash_type == ("md5-d" if is_folder else "md5")
    if is_folder:
        assert loaded.n_files == 2
        assert (loaded.path / "sub" / "b.txt").read_bytes() == b"bb"
    else:
        assert loaded.path.read_bytes() == path.read_bytes()
    assert path.exists()


def test_from_remote_reimport_returns_existing(source, target_storage):
    path = source / "file.txt"
    # memory:// has no hash, simulate a source with a checksum such as an ETag
    stats = {"size": path.stat().st_size, "hash": "source-etag", "hash_type": "md5"}
    first = ln.Artifact.from_remote(
        path, key="imports/file.txt", storage=target_storage, _stats=stats
    ).save()
    assert first.hash != "source-etag"
    second = ln.Artifact.from_remote(
        path, key="imports/file.txt", storage=target_storage, _stats=stats
    )
    assert not second._state.adding
    assert second.id == first.id
    assert second.extra_data["source"] == first.extra_data["source"]
    assert ln.Artifact.filter(storage=target_storage).count() == 1


def test_from_remote_reimport_without_source_hash(source, target_storage):
    path = source / "file.txt"
    first = ln.Artifact.from_remote(
        path, key="imports/file.txt", storage=target_storage
    ).save()
    second = ln.Artifact.from_remote(
        path, key="imports/file.txt", storage=target_storage
    ).save()
    assert second.id == first.id
    first.refresh_from_db()
    assert first.is_latest
    assert ln.Artifact.filter(storage=target_storage).count() == 1
    stored = [p for p in Path(target_storage.root).rglob("*.txt") if p.is_file()]
    assert [p.name for p in stored if p.name != "storage_uid.txt"] == [first.path.name]


def test_from_remote_local_path_matches_constructor(tmp_path, target_storage):
    path = tmp_path / "local.txt"
    path.write_text(uuid4().hex)
    direct = ln.Artifact(path, key="local/local.txt", storage=target_storage)
    imported = ln.Artifact.from_remote(
        path, key="local/local.txt", storage=target_storage
    )
    for name in ("storage_id", "key", "size", "hash", "extra_data"):
        assert getattr(imported, name) == getattr(direct, name)


def test_from_remote_rejects_dataframe():
    with pytest.raises(TypeError, match="path"):
        ln.Artifact.from_remote(pd.DataFrame({"a": [1]}), key="df.parquet")


def test_from_remote_provenance_excludes_credentials(source, target_storage):
    endpoint = "https://user:endpoint-secret@example.org:9000/api?token=query-secret"
    path = UPath(
        (source / "file.txt").as_posix(),
        endpoint_url=endpoint,
        key="access-key-secret",
        secret="secret-key-secret",  # noqa: S106
    )
    artifact = ln.Artifact.from_remote(
        path,
        key="imports/credentials.txt",
        storage=target_storage,
        extra_data={"caller": True, "source": "ignored"},
    )
    provenance = artifact.extra_data["source"]
    assert provenance["endpoint_url"] == "https://example.org:9000/api"
    assert artifact.extra_data["caller"] is True
    serialized = json.dumps(artifact.extra_data)
    for secret in ("user", "secret", "query-secret"):
        assert secret not in serialized

    https = ln.Artifact.from_remote(
        "https://user:pw@example.org/file.txt?token=abc#frag",
        description="sanitized URL",
        storage=target_storage,
        skip_check_exists=True,
        _stats={"size": 7},
    )
    assert https.extra_data["source"]["url"] == "https://example.org/file.txt"


def test_from_remote_failed_copy_leaves_no_record(source, target_storage):
    artifact = ln.Artifact.from_remote(
        source / "file.txt", key="failed/file.txt", storage=target_storage
    )
    with patch(
        "lamindb.core.storage.paths._store_from_remote",
        side_effect=OSError("interrupted"),
    ):
        with pytest.raises(OSError, match="interrupted"):
            artifact.save()
    assert not ln.Artifact.filter(uid=artifact.uid).exists()
    assert not list(Path(target_storage.root).rglob("file.txt"))


@pytest.mark.parametrize("print_progress", [False, True])
def test_from_remote_forwards_progress_and_store_kwargs(
    source, target_storage, capsys, print_progress
):
    from lamindb.core.storage import paths

    artifact = ln.Artifact.from_remote(
        source / "folder", key="imports/folder", storage=target_storage
    )
    capsys.readouterr()
    with patch.object(
        paths, "_store_from_remote", wraps=paths._store_from_remote
    ) as store:
        artifact.save(print_progress=print_progress, store_kwargs={"autocommit": True})
    assert store.call_args.kwargs == {
        "print_progress": print_progress,
        "autocommit": True,
    }
    assert ("uploading folder: 100.0%" in capsys.readouterr().out) == print_progress


def test_from_remote_copies_folder_members_in_parallel(source, target_storage):
    from concurrent.futures import ThreadPoolExecutor

    from lamindb.core.storage import paths

    for i in range(20):
        (source / "many" / f"{i}.txt").write_bytes(str(i).encode())
    artifact = ln.Artifact.from_remote(
        source / "many", key="imports/many", storage=target_storage
    )
    with patch.object(paths, "ThreadPoolExecutor", wraps=ThreadPoolExecutor) as pool:
        artifact.save(store_kwargs={"max_workers": 4})
    pool.assert_called_once_with(max_workers=4)
    assert artifact.n_files == 20
    for i in range(20):
        assert (artifact.path / f"{i}.txt").read_text() == str(i)
