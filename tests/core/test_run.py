import time

import lamindb as ln
import pytest


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
    run = ln.Run(transform).save()
    assert run.status == "scheduled"
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
    ln.Run.filter(id__in=run_ids).delete(permanent=True)
    assert ln.Run.filter(id__in=run_ids).count() == 0
    assert ln.Artifact.filter(uid=report_artifacts[0].uid).count() == 1
    transform.delete(permanent=True)

    # wait for background cleanup subprocess to delete artifacts
    time.sleep(4)
    assert ln.Artifact.filter(uid=report_artifacts[0].uid).count() == 0
