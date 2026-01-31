from unittest.mock import patch

import lamindb as ln
import lamindb_setup as ln_setup
import pytest


def test_run():
    transform = ln.Transform(key="My transform")
    with pytest.raises(ValueError) as error:
        ln.Run(transform)
    assert (
        error.exconly()
        == "ValueError: Please save transform record before creating a run"
    )
    transform.save()
    run = ln.Run(transform).save()
    assert run.status == "scheduled"
    assert run.reference is None
    assert run.reference_type is None
    run2 = ln.Run(transform, reference="test1", reference_type="test2").save()
    assert run2.reference == "test1"
    assert run2.reference_type == "test2"
    assert run.uid != run2.uid

    report_artifact = ln.Artifact("README.md", description="report of run2").save()
    run2.report = report_artifact
    environment = ln.Artifact("CONTRIBUTING.md", description="env of run2").save()
    run2.environment = environment

    run2.delete(permanent=True)

    # Run is deleted; report/env artifacts are cleaned up in background subprocess
    assert ln.Run.filter(uid=run2.uid).count() == 0
    # Clean up orphan report/env artifacts (subprocess may not have run yet)
    for art in [report_artifact, environment]:
        if ln.Artifact.filter(uid=art.uid).first() is not None:
            art.delete(permanent=True, storage=False)

    transform.delete(permanent=True)

    assert ln.Run.filter(uid=run.uid).count() == 0


def test_edge_cases():
    with pytest.raises(ValueError) as error:
        ln.Run(1, 2)
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed: transform"
    with pytest.raises(TypeError) as error:
        ln.Run()
    assert error.exconly() == "TypeError: Pass transform parameter"


def test_bulk_run_permanent_delete(tmp_path):
    """Bulk Run permanent delete uses single SQL DELETE and spawns artifact cleanup."""
    transform = ln.Transform(key="Bulk run delete transform").save()
    runs = [ln.Run(transform).save() for _ in range(3)]
    report_files = [tmp_path / f"report_{i}.txt" for i in range(3)]
    for f in report_files:
        f.write_text("report content")
    report_artifacts = [
        ln.Artifact(str(f), description=f"report {i}").save()
        for i, f in enumerate(report_files)
    ]
    for run, art in zip(runs, report_artifacts):
        run.report = art
        run.save()
    run_ids = [r.id for r in runs]
    artifact_ids = [r.report_id for r in runs]

    with patch("lamindb.models.run.subprocess.Popen") as mock_popen:
        ln.Run.filter(id__in=run_ids).delete(permanent=True)
        mock_popen.assert_called_once()
        args = mock_popen.call_args[0][0]
        assert args[args.index("--instance") + 1] == ln_setup.settings.instance.slug
        ids_str = args[args.index("--ids") + 1]
        assert {int(x) for x in ids_str.split(",")} == set(artifact_ids)

    for rid in run_ids:
        assert ln.Run.filter(id=rid).count() == 0
    # With mock, cleanup subprocess did not run; clean up orphan report artifacts
    for aid in artifact_ids:
        art = ln.Artifact.filter(id=aid).first()
        if art is not None:
            art.delete(permanent=True, storage=False)

    transform.delete(permanent=True)


def test_bulk_run_soft_delete():
    """Bulk Run soft delete sets branch_id=-1."""
    transform = ln.Transform(key="Bulk run soft delete transform").save()
    runs = [ln.Run(transform).save() for _ in range(2)]
    run_ids = [r.id for r in runs]
    ln.Run.filter(id__in=run_ids).delete(permanent=False)
    for run in ln.Run.filter(id__in=run_ids):
        assert run.branch_id == -1
    ln.Run.filter(id__in=run_ids).delete(permanent=True)
    transform.delete(permanent=True)


def test_empty_run_queryset_delete_no_subprocess():
    """Empty Run queryset delete does not spawn cleanup subprocess."""
    with patch("lamindb.models.run.subprocess.Popen") as mock_popen:
        ln.Run.filter(id=-999).delete(permanent=True)
        mock_popen.assert_not_called()
