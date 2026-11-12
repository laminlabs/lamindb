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
    from django.db import transaction
    from lamindb.core._finish import save_context_core, save_run_logs
    from lamindb_setup.core.hashing import hash_file

    with TemporaryDirectory() as directory, pytest.MonkeyPatch.context() as monkeypatch:
        tmp_path = Path(directory)
        monkeypatch.setattr("builtins.input", lambda _: "y")
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
            run = ln.Run(transform, params={"options": {"token": secret}, "safe": 42})
            run.cli_args = f"script --token {secret} --count 3"
            if bulk:
                ln.save([run])
            else:
                run.save()
            run.refresh_from_db()
            transform.refresh_from_db()
            run.params = {"options": {"token": secret}, "safe": 42}
            if bulk:
                ln.save([run])
            else:
                run.save(update_fields=["params"])
            run.refresh_from_db()
            assert run.params["safe"] == 42
            assert "--count 3" in run.cli_args

            logs = ln.setup.settings.cache_dir / f"run_logs_{run.uid}.txt"
            logs.write_text(f"logged {secret}")
            save_run_logs(run, save_run=True)
            assert secret not in run.report.path.read_text()
            assert secret in logs.read_text()
            assert hash_file(run.report.path)[1] == run.report.hash

            # Different inputs sanitize to the same content; reuse before INSERT,
            # including inside the transaction used by bulk saves.
            duplicate_path = tmp_path / "duplicate.txt"
            duplicate_path.write_text("logged ghp_" + "D" * 36)
            duplicate = ln.Artifact(
                duplicate_path,
                kind="__lamindb_run__",
                description="same log",
                run=False,
            )
            assert duplicate.pk == run.report_id
            duplicate_stage = Path(duplicate._redaction_stage.name)
            with transaction.atomic():
                ln.save([duplicate]) if bulk else duplicate.save()
            assert not duplicate_stage.exists()

            script = tmp_path / "script.py"
            script.write_text('print("ordinary_value")')
            environment = (
                ln.setup.settings.cache_dir / "environments" / f"run_{run.uid}"
            )
            environment.mkdir(parents=True)
            for name in ("uv.lock", "pixi.lock", "run_env_pip.txt"):
                (environment / name).write_text(
                    f"pkg @ https://user:{secret}@example.org/pkg"
                )
            save_context_core(
                run=run, transform=transform, filepath=script, from_cli=True
            )
            run.refresh_from_db()
            assert secret in (environment / "uv.lock").read_text()
            # Replacement uses the same protection and hashes the redacted bytes.
            replacement = tmp_path / "replacement.txt"
            replacement.write_text(f"new log {secret}")
            run.report.replace(replacement)
            stage = Path(run.report._redaction_stage.name)
            ln.save([run.report]) if bulk else run.report.save()
            assert not stage.exists()
            assert not hasattr(run.report, "_redaction_stage")
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
            ln.Artifact(source, kind="__lamindb_run__", description="report", run=False)
        report = ln.Artifact(
            source,
            kind="__lamindb_run__",
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
def test_run_attachment_rejects_changed_snapshot(tmp_path, bulk):
    from lamindb.core._secret_redaction import _UnsafeRunArtifact

    source = tmp_path / "log.txt"
    source.write_text("original log")
    artifact = ln.Artifact(
        source, kind="__lamindb_run__", description="run log", run=False
    )
    artifact._local_filepath.write_text('api_key = "CANARY_CHANGED"')
    with pytest.raises(_UnsafeRunArtifact, match="changed after preparation"):
        ln.save([artifact]) if bulk else artifact.save()
    assert not ln.Artifact.filter(uid=artifact.uid).exists()
    assert source.read_text() == "original log"
    artifact._redaction_stage.cleanup()
