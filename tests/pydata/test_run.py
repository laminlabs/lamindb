import time
from pathlib import Path
from string import ascii_letters, digits
from tempfile import TemporaryDirectory

import lamindb as ln
import pytest
from hypothesis import given, settings
from hypothesis import strategies as st


@pytest.mark.parametrize("suffix", ["", ".html", ".txt"])
@pytest.mark.parametrize("is_ipython", [False, True])
def test_describe_run_report(tmp_path, monkeypatch, suffix, is_ipython):
    from IPython.display import HTML
    from lamindb.core import loaders

    displayed = []
    monkeypatch.setattr(loaders, "is_run_from_ipython", is_ipython)
    monkeypatch.setattr("IPython.display.display", displayed.append)
    report_path = tmp_path / f"report{suffix}"
    content = "run report content"
    if suffix == "":
        report_path.mkdir()
        html_path = report_path / "index.html"
        (report_path / "trace.jsonl").write_text('{"role": "user"}\n')
    else:
        html_path = report_path
    html_path.write_text(
        content if suffix == ".txt" else f"<html><body>{content}</body></html>"
    )
    artifact = ln.Artifact(report_path, description="run report", run=False).save()
    transform = ln.Transform(key="describe_run_report.py").save()
    run = ln.Run(transform, report=artifact).save()
    try:
        description = run.describe(return_str=True)
        if suffix == ".txt":
            assert content in description
            assert displayed == []
        elif is_ipython:
            assert len(displayed) == 1
            assert isinstance(displayed[0], HTML)
            assert content in displayed[0].data
        else:
            assert f"report: {artifact.uid[:7]}" in description
            assert displayed == []
    finally:
        run.report = None
        run.save()
        run.delete(permanent=True)
        artifact.delete(permanent=True)
        transform.delete(permanent=True)


def _wait_for_artifact_cleanup(*artifact_uids: str, timeout: float = 10.0) -> None:
    """Poll until report/env artifacts are deleted by the background cleanup subprocess.

    That subprocess cold-imports lamindb, which on CI (esp. Python 3.12 with --cov)
    can take longer than a fixed sleep.
    """
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        if not ln.Artifact.filter(uid__in=artifact_uids).exists():
            return
        time.sleep(0.5)


def test_run():
    with pytest.raises(ValueError) as error:
        ln.Run(1, 2)
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed: transform"
    with pytest.raises(TypeError) as error:
        ln.Run()
    assert error.exconly() == "TypeError: Pass transform parameter"
    transform = ln.Transform(key="my_transform")
    with pytest.raises(ValueError) as error:
        ln.Run(transform)
    assert (
        error.exconly()
        == "ValueError: Please save transform record before creating a run"
    )
    transform.save()
    with pytest.raises(ValueError) as error:
        ln.Run(transform, status="invalid")
    assert (
        error.exconly()
        == "ValueError: status must be 'scheduled' or 'started', but you passed: 'invalid'."
    )
    run = ln.Run(transform).save()
    assert run.status == "scheduled"
    started_run = ln.Run(transform, status="started").save()
    assert started_run.status == "started"
    assert run.reference is None
    assert run.reference_type is None
    run2 = ln.Run(transform, reference="test1", reference_type="test2").save()
    assert run2.reference == "test1"
    assert run2.reference_type == "test2"
    assert run.uid != run2.uid
    run.delete(permanent=True)

    report_artifact = ln.Artifact(
        "README.md", kind="__lamindb_run__", description="report of run2"
    ).save()
    run2.report = report_artifact
    environment = ln.Artifact(
        "CONTRIBUTING.md", kind="__lamindb_run__", description="requirements.txt"
    ).save()
    run2.environment = environment
    run2.save()

    # report/env artifacts will be cleaned up in background subprocess
    run2.delete(permanent=True)
    assert ln.Run.filter(uid=run2.uid).count() == 0
    # report/env are still present in the database
    assert ln.Artifact.filter(uid=report_artifact.uid).count() == 1
    assert ln.Artifact.filter(uid=environment.uid).count() == 1

    transform.delete(permanent=True)
    assert ln.Run.filter(uid=run.uid).count() == 0

    # wait for background cleanup subprocess to delete artifacts
    time.sleep(4)
    _wait_for_artifact_cleanup(report_artifact.uid, environment.uid)
    assert ln.Artifact.filter(uid=report_artifact.uid).count() == 0
    assert ln.Artifact.filter(uid=environment.uid).count() == 0


def test_bulk_permanent_run_delete(tmp_path):
    transform = ln.Transform(key="Bulk run delete transform").save()
    n_runs = 2
    report_files = [tmp_path / f"report_{i}.txt" for i in range(n_runs)]
    for i, path in enumerate(report_files):
        path.write_text(f"content {i}")
    report_artifacts = [
        ln.Artifact(path, kind="__lamindb_run__", description=f"report {i}").save()
        for i, path in enumerate(report_files)
    ]
    runs = [ln.Run(transform, report=af).save() for af in report_artifacts]
    run_ids = [r.id for r in runs]
    ln.settings.verbosity = "debug"
    ln.Run.filter(id__in=run_ids).order_by("created_at").delete(permanent=True)
    assert ln.Run.filter(id__in=run_ids).count() == 0
    assert ln.Artifact.filter(uid=report_artifacts[0].uid).count() == 1
    transform.delete(permanent=True)

    # wait for background cleanup subprocess to delete artifacts
    time.sleep(4)
    _wait_for_artifact_cleanup(report_artifacts[0].uid)
    assert ln.Artifact.filter(uid=report_artifacts[0].uid).count() == 0
    clean_up_logs = ln.setup.settings.cache_dir / f"run_cleanup_logs_{runs[0].uid}.txt"
    assert f"deleted artifact {report_artifacts[0].id}" in clean_up_logs.read_text()


@pytest.mark.parametrize("bulk", [False, True])
@settings(max_examples=20, deadline=None)
@given(
    secret=st.text(ascii_letters + digits, min_size=36, max_size=36).map(
        "ghp_{}".format
    )
)
def test_stored_run_secrets(bulk, secret):
    from datetime import UTC, datetime
    from types import SimpleNamespace

    from django.db import transaction
    from lamindb.core._finish import save_context_core, save_run_logs
    from lamindb.core._secret_redaction import redact
    from lamindb.core._track_environment import _track_pip_freeze
    from lamindb_setup.core.hashing import hash_file

    with TemporaryDirectory() as directory, pytest.MonkeyPatch.context() as monkeypatch:
        tmp_path = Path(directory)
        monkeypatch.setattr("builtins.input", lambda _: "y")
        monkeypatch.setattr(ln.settings.creation, "artifact_skip_size_hash", bulk)
        before = set(
            ln.Artifact.filter(kind="__lamindb_run__").values_list("id", flat=True)
        )
        transform = ln.Transform(
            key=f"redaction-{bulk}.py", source_code=f'x = "{secret}"'
        )
        run = None
        try:
            if bulk:
                ln.save([transform])
            else:
                transform.save()
            params = {
                "options": {"token": secret},
                "labels": {secret: "primary"},
                secret: "top level",
                "safe": 42,
            }
            run = ln.Run(transform, params=params)
            run.cli_args = f"script --token {secret} --count 3"
            if bulk:
                ln.save([run])
            else:
                run.save()
            run.refresh_from_db()
            transform.refresh_from_db()
            run.params = params
            if bulk:
                ln.save([run])
            else:
                run.save(update_fields=["params"])
            run.refresh_from_db()
            assert run.params["safe"] == 42
            assert run.params["labels"] == {"***REDACTED***": "primary"}
            assert run.params["***REDACTED***"] == "top level"
            assert "--count 3" in run.cli_args

            logs = ln.setup.settings.cache_dir / f"run_logs_{run.uid}.txt"
            logs.write_text(f"logged {secret}")
            run._status_code = status_code = 2 if bulk else 1
            run.finished_at = finished_at = datetime.now(UTC)
            save_run_logs(run, save_run=True)
            run.refresh_from_db()
            assert run._status_code == status_code
            assert run.finished_at == finished_at
            assert run.report is None  # Existing unsafe logs are skipped, not copied.
            assert secret in logs.read_text()
            logs.write_text("ordinary log")
            save_run_logs(run, save_run=True)
            assert run.report.path.read_text() == "ordinary log"
            assert hash_file(run.report.path)[1] == run.report.hash

            # Already-safe files reuse content before INSERT,
            # including inside the transaction used by bulk saves.
            duplicate_path = tmp_path / "duplicate.txt"
            duplicate_path.write_text("ordinary log")
            duplicate = ln.Artifact(
                duplicate_path,
                kind="__lamindb_run__",
                _check_run=True,
                description="same log",
                run=False,
            )
            assert duplicate.pk == run.report_id
            with transaction.atomic():
                ln.save([duplicate]) if bulk else duplicate.save()

            script = tmp_path / "script.py"
            script.write_text('print("ordinary_value")')
            environment = (
                ln.setup.settings.cache_dir / "environments" / f"run_{run.uid}"
            )
            monkeypatch.setattr(
                "lamindb.core._track_environment.subprocess.run",
                lambda *a, **kw: SimpleNamespace(
                    returncode=0, stdout=f"pkg @ https://user:{secret}@example.org/pkg"
                ),
            )
            assert _track_pip_freeze(environment) is None
            save_context_core(
                run=run, transform=transform, filepath=script, from_cli=True
            )
            run.refresh_from_db()
            assert secret not in (environment / "run_env_pip.txt").read_text()
            # Generated replacements are redacted before writing and hashing.
            replacement = tmp_path / "replacement.txt"
            replacement.write_text(redact(f"new log {secret}"))
            run.report.replace(replacement, _check_run=True)
            ln.save([run.report]) if bulk else run.report.save()
            assert hash_file(run.report.path)[1] == run.report.hash
            # Scan all persisted columns and every new run attachment, not a field allowlist.
            for record in (run, transform):
                stored = type(record).objects.filter(pk=record.pk).values().get()
                assert secret not in str(stored)
            artifacts = list(
                ln.Artifact.filter(kind="__lamindb_run__").exclude(id__in=before)
            )
            assert {run.report_id, run.environment_id} <= {a.id for a in artifacts}
            for artifact in artifacts:
                path = artifact.path
                for file in path.rglob("*") if path.is_dir() else [path]:
                    if file.is_file():
                        assert secret.encode() not in file.read_bytes()
        finally:
            artifacts = list(
                ln.Artifact.filter(kind="__lamindb_run__").exclude(id__in=before)
            )
            if run is not None and not run._state.adding:
                run.delete(permanent=True)
            if not transform._state.adding:
                transform.delete(permanent=True)
            for artifact in artifacts:
                artifact.delete(permanent=True, storage=True)


@pytest.mark.parametrize("content", ["clean log", ""])
def test_run_attachment_dedup_preserves_kind_and_storage(tmp_path, content):
    from lamindb.core._secret_redaction import _UnsafeRunArtifact

    source = tmp_path / "source.txt"
    source.write_text(content)
    ordinary = ln.Artifact(source, description="ordinary data", run=False).save()
    storage = ln.Storage(root=(tmp_path / "other-storage").as_posix()).save()
    report = None
    try:
        with pytest.raises(_UnsafeRunArtifact, match="incompatible artifact"):
            ln.Artifact(
                source,
                kind="__lamindb_run__",
                _check_run=True,
                description="report",
                run=False,
            )
        report = ln.Artifact(
            source,
            kind="__lamindb_run__",
            _check_run=True,
            description="report",
            storage=storage,
            run=False,
        ).save()
        assert report.pk != ordinary.pk and report.storage_id == storage.id
        ordinary.refresh_from_db()
        assert ordinary.kind != "__lamindb_run__"
        duplicate = ln.Artifact(
            source,
            kind="__lamindb_run__",
            _check_run=True,
            description="report",
            storage=storage,
            run=False,
        ).save()
        assert duplicate.pk == report.pk
    finally:
        if report is not None:
            report.delete(permanent=True, storage=True)
        ordinary.delete(permanent=True, storage=True)
        storage.delete()


@pytest.mark.parametrize("bulk", [False, True])
def test_run_attachment_rejects_changed_file(tmp_path, bulk):
    from lamindb.core._secret_redaction import _UnsafeRunArtifact

    source = tmp_path / "log.txt"
    source.write_text("original log")
    artifact = ln.Artifact(
        source,
        kind="__lamindb_run__",
        _check_run=True,
        description="run log",
        run=False,
    )
    artifact._local_filepath.write_text("changed clean log")
    with pytest.raises(_UnsafeRunArtifact, match="changed after preparation"):
        ln.save([artifact]) if bulk else artifact.save()
    assert not ln.Artifact.filter(uid=artifact.uid).exists()
    assert artifact._local_filepath == source
    assert source.read_text() == "changed clean log"


@pytest.mark.parametrize(
    "match", ["run", "ordinary", "other_storage", "other_space", "incomplete"]
)
def test_run_attachment_replace_duplicate(tmp_path, match):
    from lamindb.core._secret_redaction import _UnsafeRunArtifact

    source = tmp_path / "old.txt"
    source.write_text("old report")
    replacement = tmp_path / "new.txt"
    replacement.write_text('new report: api_key = "***REDACTED***"')
    other_storage = ln.Storage(root=(tmp_path / "other").as_posix()).save()
    report = ln.Artifact(
        source,
        kind="__lamindb_run__",
        _check_run=True,
        description="report",
        run=False,
        storage=other_storage,
    ).save()
    # Ordinary artifacts need already-sanitized content to collide.
    duplicate_path = tmp_path / "duplicate.txt"
    duplicate_path.write_text('new report: api_key = "***REDACTED***"')
    duplicate = ln.Artifact(
        duplicate_path,
        kind=None if match == "ordinary" else "__lamindb_run__",
        description="duplicate",
        run=False,
        storage=None if match == "other_storage" else other_storage,
    ).save()
    space = None
    try:
        if match == "other_space":
            space = ln.Space(name="replacement-space").save()
            ln.Artifact.objects.filter(pk=duplicate.pk).update(space=space)
        elif match == "incomplete":
            duplicate._storage_ongoing = True
            ln.Artifact.objects.filter(pk=duplicate.pk).update(_aux=duplicate._aux)
        if match in {"ordinary", "other_space"}:
            with pytest.raises(_UnsafeRunArtifact, match="incompatible artifact"):
                report.replace(replacement, run=False, _check_run=True)
            assert report.path.read_text() == "old report"
        else:
            existing = report.replace(replacement, run=False, _check_run=True)
            if match == "other_storage":
                assert existing is None
                report.save()
                updated = report
            else:
                assert existing.pk == duplicate.pk
                existing.save()
                updated = existing
            assert updated.storage_id == other_storage.id
            assert updated.path.read_text() == duplicate_path.read_text()
        assert replacement.read_text() == duplicate_path.read_text()
    finally:
        report.delete(permanent=True, storage=True)
        duplicate.delete(permanent=True, storage=True)
        if space is not None:
            space.delete()
        other_storage.delete()


def test_ordinary_artifact_still_skips_hash(tmp_path, monkeypatch):
    monkeypatch.setattr(ln.settings.creation, "artifact_skip_size_hash", True)
    path = tmp_path / "ordinary.txt"
    path.write_text("ordinary data")
    artifact = ln.Artifact(path, description="ordinary", run=False).save()
    try:
        artifact.refresh_from_db()
        assert artifact.size is None and artifact.hash is None
    finally:
        artifact.delete(permanent=True, storage=True)


@pytest.mark.parametrize("set_kind_after_init", [False, True])
def test_internal_binary_artifact_is_not_a_run_attachment(
    tmp_path, set_kind_after_init
):
    import pandas as pd

    path = tmp_path / "ontology.parquet"
    frame = pd.DataFrame({"name": ["T cell"], "ontology_id": ["CL:0000084"]})
    frame.to_parquet(path)
    original = path.read_bytes()
    artifact = (
        ln.Artifact.from_dataframe(path, description="ontology", run=False)
        if set_kind_after_init
        else ln.Artifact(
            path, kind="__lamindb_run__", description="ontology", run=False
        )
    )
    # Bionty marks reference data as internal after constructing the artifact.
    artifact.kind = "__lamindb_run__"
    artifact.save()
    try:
        assert artifact.path.read_bytes() == original
        pd.testing.assert_frame_equal(pd.read_parquet(artifact.path), frame)
    finally:
        artifact.delete(permanent=True, storage=True)


@pytest.mark.parametrize("existing_target", [False, True])
@pytest.mark.parametrize("custom_space", [False, True])
def test_run_environment_reuse_scope(
    tmp_path, monkeypatch, existing_target, custom_space
):
    from lamindb.core._finish import save_context_core

    space = ln.Space(name="environment-reuse-space").save()
    other_storage = ln.Storage(root=(tmp_path / "other").as_posix()).save()
    # Set up local space metadata without registering a space with the Hub.
    ln.Storage.objects.filter(pk=other_storage.pk).update(space=space)
    other_storage.refresh_from_db()
    default_storage = ln.setup.settings.storage.record
    target, foreign = (
        (other_storage, default_storage)
        if custom_space
        else (default_storage, other_storage)
    )
    transform = ln.Transform(key="environment-reuse.py").save()
    run = ln.Run(transform).save()
    artifacts = []
    try:
        seed = tmp_path / "seed.txt"
        seed.write_text("pkg @ https://user:***REDACTED***@example.org/pkg")
        for storage in [foreign, target] if existing_target else [foreign]:
            artifacts.append(
                ln.Artifact(
                    seed,
                    storage=storage,
                    kind="__lamindb_run__",
                    _check_run=True,
                    description="environment",
                    run=False,
                ).save()
            )
        monkeypatch.setattr(ln.context, "_space", space if custom_space else None)
        environment = ln.setup.settings.cache_dir / "environments" / f"run_{run.uid}"
        environment.mkdir(parents=True)
        source = environment / "run_env_pip.txt"
        source.write_text("pkg @ https://user:***REDACTED***@example.org/pkg")
        script = tmp_path / "script.py"
        script.write_text('print("ordinary_value")')
        save_context_core(run=run, transform=transform, filepath=script, from_cli=True)
        run.refresh_from_db()
        artifact = run.environment
        if not existing_target:
            artifacts.append(artifact)
        else:
            assert artifact.pk == artifacts[-1].pk
        assert artifact.storage_id == target.id
        assert artifact.space_id == target.space_id
        assert artifact.pk != artifacts[0].pk
        assert "CANARY" not in artifact.path.read_text()
        assert source.read_text() == artifact.path.read_text()
    finally:
        run.environment = None
        run.save()
        run.delete(permanent=True)
        transform.delete(permanent=True)
        for artifact in artifacts:
            artifact.delete(permanent=True, storage=True)
        other_storage.delete()
        space.delete()
