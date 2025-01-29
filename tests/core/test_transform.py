from pathlib import Path

import lamindb as ln
import pytest


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
    assert transform.version == "1"
    assert len(transform.uid) == ln.Transform._len_full_uid == 16
    assert len(transform.stem_uid) == ln.Transform._len_stem_uid == 12

    transform.save()

    # try to reload the same transform with the same uid
    transform_reload = ln.Transform(uid=transform.uid, name="My transform updated name")
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
    assert transform_r2.version is None
    assert transform_r2.is_latest
    assert transform.is_latest
    transform_r2.save()
    assert not transform.is_latest

    # create new transform from newly versioned transform
    transform_r3 = ln.Transform(
        description="My transform", revises=transform_r2, version="2"
    )
    assert transform_r3.stem_uid == transform.stem_uid
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
        ln.Transform(revises=ln.ULabel(name="x"))
    assert error.exconly().startswith(
        "TypeError: `revises` has to be of type `Transform`"
    )

    # wrong kwargs
    with pytest.raises(ValueError) as error:
        ln.Transform(x=1)
    assert (
        error.exconly() == "ValueError: Only key, description, version, type, revises,"
        " reference, reference_type can be passed, but you passed: {'x': 1}"
    )

    # test that reference transform cannot be deleted
    transform_r2.delete()
    transform.delete()

    # unversioned transform
    transform = ln.Transform(key="My transform")
    assert transform.version is None

    # what happens if we don't save the old transform?
    # add a test for it!
    transform.save()

    # create new transform from old transform
    new_transform = ln.Transform(description="My new transform", revises=transform)
    assert transform.version is None
    assert new_transform.stem_uid == transform.stem_uid
    assert new_transform.uid.endswith("0001")
    assert new_transform.version is None

    transform.delete()


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
    # now delete everything
    transform.delete()
    assert not report_path.exists()
    assert not environment_path.exists()
    assert len(ln.Artifact.filter(id__in=[report.id, environment.id]).all()) == 0
    assert len(ln.Run.filter(id=run.id).all()) == 0
