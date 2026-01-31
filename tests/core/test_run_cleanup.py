"""Tests for Run bulk delete and _run_cleanup subprocess."""

import subprocess
import sys

import lamindb as ln
import lamindb_setup as ln_setup


def test_run_cleanup_module_nonexistent_id_no_crash():
    """Running _run_cleanup with non-existent artifact ID does not crash."""
    instance = ln_setup.settings.instance.slug
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "lamindb.models._run_cleanup",
            "--instance",
            instance,
            "--ids",
            "999999",
        ],
        capture_output=True,
        timeout=10,
    )
    assert result.returncode == 0


def test_run_cleanup_module_orphan_artifact_deleted(tmp_path):
    """Running _run_cleanup with orphan artifact ID deletes the artifact."""
    # Create an artifact that is not linked to any run (orphan)
    f = tmp_path / "orphan_cleanup_test.txt"
    f.write_text("orphan")
    art = ln.Artifact(str(f), description="orphan for cleanup").save()
    art_id = art.id
    assert ln.Artifact.filter(id=art_id).first() is not None

    instance = ln_setup.settings.instance.slug
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "lamindb.models._run_cleanup",
            "--instance",
            instance,
            "--ids",
            str(art_id),
        ],
        capture_output=True,
        timeout=10,
    )
    assert result.returncode == 0
    assert ln.Artifact.filter(id=art_id).first() is None
