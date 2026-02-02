from pathlib import Path
from unittest.mock import patch

import lamindb as ln
import pytest


def test_transform_recovery_based_on_hash():
    transform1 = ln.Transform(key="my-transform", source_code="1").save()
    transform2 = ln.Transform(key="my-transform", source_code="1")
    assert transform1 == transform2
    transform1.delete()
    transform2 = ln.Transform(key="my-transform", source_code="1")
    assert transform1 != transform2
    transform1.delete(permanent=True)


def test_transform_recovery_based_on_key():
    transform1 = ln.Transform(key="my-transform").save()
    transform2 = ln.Transform(key="my-transform")
    assert transform1 == transform2
    transform1.delete()
    transform2 = ln.Transform(key="my-transform")
    assert transform1 != transform2
    transform1.delete(permanent=True)


def test_revise_transforms():
    # attempt to create a transform with an invalid version
    with pytest.raises(ValueError) as error:
        transform = ln.Transform(key="My transform", version=0)
        assert (
            error.exconly()
            == "ValueError: `version` parameter must be `None` or `str`, e.g., '0.1', '1',"
            " '2', etc."
        )

    # create a versioned transform
    transform = ln.Transform(key="My transform", version="1")
    assert transform.version_tag == "1"
    assert transform.version == "1"
    assert len(transform.uid) == ln.Transform._len_full_uid == 16
    assert len(transform.stem_uid) == ln.Transform._len_stem_uid == 12

    transform.save()

    # try to reload the same transform with the same uid
    transform_reload = ln.Transform(uid=transform.uid, key="My transform updated name")
    assert transform_reload.id == transform.id
    assert transform_reload.key == "My transform"  # unchanged, prints logging
    transform_reload = ln.Transform(
        uid=transform.uid, description="My transform updated name"
    )
    assert transform_reload.id == transform.id
    assert (
        transform_reload.description == "My transform updated name"
    )  # unchanged, prints logging

    # create new transform from old transform
    transform_r2 = ln.Transform(description="My 2nd transform", revises=transform)
    assert transform_r2.uid != transform.uid
    assert transform_r2.uid.endswith("0001")
    transform_r2 = ln.Transform(description="My 2nd transform", revises=transform)
    assert transform_r2.uid != transform.uid
    assert transform_r2.uid.endswith("0001")
    assert transform_r2.stem_uid == transform.stem_uid
    assert transform_r2.version_tag is None
    assert (
        transform_r2.version == transform_r2.uid[-4:]
    )  # version falls back to uid suffix
    assert transform_r2.is_latest
    assert transform.is_latest
    transform_r2.save()
    assert not transform.is_latest

    # create new transform from newly versioned transform
    transform_r3 = ln.Transform(
        description="My transform", revises=transform_r2, version="2"
    )
    assert transform_r3.stem_uid == transform.stem_uid
    assert transform_r3.version_tag == "2"
    assert transform_r3.version == "2"

    # default description
    transform_r3 = ln.Transform(revises=transform_r2)
    assert transform_r3.description == transform_r2.description

    # revise by matching on `key`
    key = "my-notebook.ipynb"
    transform_r2.key = key
    transform_r2.save()
    assert transform_r2.is_latest
    transform_r3 = ln.Transform(description="My transform", key=key, version="2")
    assert transform_r3.uid[:-4] == transform_r2.uid[:-4]
    assert transform_r3.uid.endswith("0001")
    # this only fires if source code was actually saved
    transform_r2.source_code = "something"
    transform_r2.save()
    transform_r3 = ln.Transform(description="My transform", key=key, version="2")
    assert transform_r3.uid[:-4] == transform_r2.uid[:-4]
    assert transform_r3.uid.endswith("0002")
    assert transform_r3.stem_uid == transform_r2.stem_uid
    assert transform_r3.key == key
    assert transform_r3.version_tag == "2"
    assert transform_r3.version == "2"
    assert transform_r3.is_latest
    # because the new transform isn't yet saved, the old transform still has
    # is_latest = True
    assert transform_r2.is_latest
    assert transform_r3._revises is not None
    transform_r3.save()
    # now r2 is no longer the latest version, but need to re-fresh from db
    transform_r2 = ln.Transform.get(transform_r2.uid)
    assert not transform_r2.is_latest

    # wrong transform type
    with pytest.raises(TypeError) as error:
        ln.Transform(revises=ln.Record(name="x"))
    assert error.exconly().startswith(
        "TypeError: `revises` has to be of type `Transform`"
    )

    # wrong kwargs
    with pytest.raises(ValueError) as error:
        ln.Transform(x=1)
        assert (
            error.exconly()
            == "ValueError: Only key, description, version_tag, type, revises,"
            " reference, reference_type can be passed, but you passed: {'x': 1}"
        )

    # test that reference transform cannot be deleted
    transform_r2.delete()
    transform.delete()

    # unversioned transform
    transform = ln.Transform(key="My transform")
    assert transform.version_tag is None
    assert transform.version == transform.uid[-4:]  # version falls back to uid suffix

    # what happens if we don't save the old transform?
    # add a test for it!
    transform.save()

    # create new transform from old transform
    new_transform = ln.Transform(description="My new transform", revises=transform)
    assert transform.version_tag is None
    assert transform.version == transform.uid[-4:]  # version falls back to uid suffix
    assert new_transform.stem_uid == transform.stem_uid
    assert new_transform.uid.endswith("0001")
    assert new_transform.version_tag is None
    assert (
        new_transform.version == new_transform.uid[-4:]
    )  # version falls back to uid suffix

    transform.delete(permanent=True)


def test_delete():
    # prepare the creation of a transform with its artifacts
    transform = ln.Transform(key="My transform").save()
    run = ln.Run(transform)
    report_path = Path("report.html")
    with open(report_path, "w") as f:
        f.write("a")
    environment_path = Path("environment.txt")
    with open(environment_path, "w") as f:
        f.write("c")
    report = ln.Artifact(report_path, description=f"Report of {run.uid}").save()
    report_path.unlink()
    report_path = report.path
    environment = ln.Artifact(environment_path, description="requirements.txt").save()
    environment_path.unlink()
    environment_path = environment.path
    transform.save()
    run.report = report
    run.environment = environment
    run.save()
    assert report_path.exists()
    assert environment_path.exists()
    # now delete everything (run artifacts are cleaned up in background subprocess)
    transform.delete(permanent=True)
    assert len(ln.Run.filter(id=run.id)) == 0
    # Clean up orphan report/env artifacts if subprocess has not run yet
    for art in [report, environment]:
        a = ln.Artifact.filter(id=art.id).first()
        if a is not None:
            a.delete(permanent=True, storage=True)
    assert not report_path.exists()
    assert not environment_path.exists()
    assert len(ln.Artifact.filter(id__in=[report.id, environment.id])) == 0


# see test_composite_component in test_schema.py
def test_successor_predecessor():
    predecessor = ln.Transform(key="predecessor").save()
    successor1 = ln.Transform(key="successor1").save()
    successor2 = ln.Transform(key="successor2").save()
    predecessor.successors.add(
        successor1, successor2, through_defaults={"config": {"param": 42}}
    )

    assert len(predecessor.successors.all()) == 2
    assert predecessor.links_successor.count() == 2
    assert predecessor.links_successor.first().config == {"param": 42}
    assert predecessor.links_successor.first().predecessor == predecessor
    assert predecessor.predecessors.count() == 0
    assert predecessor.links_predecessor.count() == 0

    ln.models.transform.TransformTransform.filter(predecessor=predecessor).delete(
        permanent=True
    )

    link = ln.models.transform.TransformTransform(
        predecessor=predecessor, successor=successor1, config={"param": 42}
    ).save()
    assert link in predecessor.links_successor.all()
    assert link in successor1.links_predecessor.all()
    assert link.config == {"param": 42}

    predecessor.delete(permanent=True)
    successor1.delete(permanent=True)
    successor2.delete(permanent=True)

    assert ln.models.transform.TransformTransform.filter().count() == 0


def test_bulk_transform_permanent_delete(tmp_path):
    """Bulk Transform permanent delete deletes TransformProject, runs (and artifacts), then transforms."""
    transform = ln.Transform(key="Bulk transform delete").save()
    runs = [ln.Run(transform).save() for _ in range(2)]
    report_files = [tmp_path / f"bulk_report_{i}.txt" for i in range(2)]
    for f in report_files:
        f.write_text("report content")
    report_artifacts = [
        ln.Artifact(str(f), description=f"report {i}").save()
        for i, f in enumerate(report_files)
    ]
    for run, art in zip(runs, report_artifacts):
        run.report = art
        run.save()
    transform_id = transform.id
    run_ids = [r.id for r in runs]
    artifact_ids = [r.report_id for r in runs]

    with patch("lamindb.models.run.subprocess.Popen") as mock_popen:
        ln.Transform.filter(id=transform_id).delete(permanent=True)
        mock_popen.assert_called_once()
        args = mock_popen.call_args[0][0]
        ids_str = args[args.index("--ids") + 1]
        assert {int(x) for x in ids_str.split(",")} == set(artifact_ids)

    assert ln.Transform.filter(id=transform_id).count() == 0
    for rid in run_ids:
        assert ln.Run.filter(id=rid).count() == 0
    # With mock, cleanup subprocess did not run; clean up orphan report artifacts
    for aid in artifact_ids:
        art = ln.Artifact.filter(id=aid).first()
        if art is not None:
            art.delete(permanent=True, storage=False)


def test_single_transform_permanent_delete_delegates_to_queryset(tmp_path):
    """Single Transform permanent delete delegates to QuerySet and removes runs and artifacts."""
    transform = ln.Transform(key="Single transform delete").save()
    run = ln.Run(transform).save()
    report_file = tmp_path / "single_report.txt"
    report_file.write_text("report")
    report = ln.Artifact(str(report_file), description="report").save()
    run.report = report
    run.save()
    transform_id = transform.id
    run_id = run.id
    artifact_id = report.id

    with patch("lamindb.models.run.subprocess.Popen") as mock_popen:
        transform.delete(permanent=True)
        mock_popen.assert_called_once()
        args = mock_popen.call_args[0][0]
        ids_str = args[args.index("--ids") + 1]
        assert artifact_id in {int(x) for x in ids_str.split(",")}

    assert ln.Transform.filter(id=transform_id).count() == 0
    assert ln.Run.filter(id=run_id).count() == 0
    # With mock, cleanup subprocess did not run; clean up orphan report artifact
    art = ln.Artifact.filter(id=artifact_id).first()
    if art is not None:
        art.delete(permanent=True, storage=False)


def test_bulk_transform_soft_delete():
    """Bulk Transform soft delete sets branch_id=-1."""
    transform = ln.Transform(key="Bulk transform soft delete").save()
    ln.Run(transform).save()
    transform_id = transform.id
    ln.Transform.filter(id=transform_id).delete(permanent=False)
    t = ln.Transform.filter(id=transform_id).one()
    assert t.branch_id == -1
    ln.Transform.filter(id=transform_id).delete(permanent=True)


def test_bulk_transform_permanent_delete_promotes_previous_version():
    """Bulk permanent delete of latest in a version family promotes the previous version."""
    v1 = ln.Transform(key="Bulk permanent delete version family").save()
    v2 = ln.Transform(revises=v1, key="Bulk permanent delete version family").save()
    assert v2.is_latest
    stem_uid = v1.stem_uid

    ln.Transform.filter(id=v2.id).delete(permanent=True)

    assert ln.Transform.filter(id=v2.id).count() == 0
    v1_after = ln.Transform.filter(uid__startswith=stem_uid).one()
    assert v1_after.pk == v1.pk
    assert v1_after.is_latest
    v1.delete(permanent=True)


def test_bulk_transform_soft_delete_promotes_previous_version():
    """Bulk soft delete of latest in a version family promotes the previous version."""
    v1 = ln.Transform(key="Bulk soft delete version family").save()
    v2 = ln.Transform(revises=v1, key="Bulk soft delete version family").save()
    assert v2.is_latest
    v2_id = v2.id
    stem_uid = v1.stem_uid

    ln.Transform.filter(id=v2_id).delete(permanent=False)

    v2_after = ln.Transform.filter(id=v2_id).one()
    assert v2_after.branch_id == -1
    assert not v2_after.is_latest
    v1.refresh_from_db()
    assert v1.is_latest
    assert ln.Transform.filter(uid__startswith=stem_uid).get(is_latest=True) == v1
    # Clean up
    v2_after.delete(permanent=True)
    v1.delete(permanent=True)
