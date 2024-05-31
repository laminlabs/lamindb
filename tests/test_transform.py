from pathlib import Path

import lamindb as ln
import pytest
from django.db.models.deletion import ProtectedError


def test_is_new_version_of_versioned_transform():
    # attempt to create a transform with an invalid version
    with pytest.raises(ValueError) as error:
        transform = ln.Transform(name="My transform", version=0)
    assert (
        error.exconly()
        == "ValueError: `version` parameter must be `None` or `str`, e.g., '0.1', '1',"
        " '2', etc."
    )

    # create a versioned transform
    transform = ln.Transform(name="My transform", version="1")
    assert transform.version == "1"
    assert len(transform.uid) == ln.Transform._len_full_uid == 16
    assert len(transform.stem_uid) == ln.Transform._len_stem_uid == 12

    transform.save()

    # create new transform from old transform
    transform_v2 = ln.Transform(name="My 2nd transform", is_new_version_of=transform)
    assert transform.version == "1"
    assert transform_v2.uid != transform.uid
    assert transform_v2.stem_uid == transform.stem_uid
    assert transform_v2.version == "2"

    transform_v2.save()

    # create new transform from newly versioned transform
    transform_v3 = ln.Transform(name="My transform", is_new_version_of=transform_v2)
    assert transform_v3.stem_uid == transform.stem_uid
    assert transform_v3.version == "3"

    # default name
    transform_v3 = ln.Transform(is_new_version_of=transform_v2)
    assert transform_v3.name == transform_v2.name

    # wrong transform type
    with pytest.raises(TypeError) as error:
        ln.Transform(is_new_version_of=ln.ULabel(name="x"))
    assert error.exconly().startswith("TypeError: is_new_version_of has to be of type")

    # wrong kwargs
    with pytest.raises(ValueError) as error:
        ln.Transform(x=1)
    assert (
        error.exconly()
        == "ValueError: Only name, key, version, type, is_new_version_of,"
        " reference, reference_type can be passed, but you passed: {'x': 1}"
    )

    # test that reference transform cannot be deleted
    transform_v2.delete()
    transform.delete()


def test_is_new_version_of_unversioned_transform():
    # unversioned transform
    transform = ln.Transform(name="My transform")
    assert transform.version is None

    # what happens if we don't save the old transform?
    # add a test for it!
    transform.save()

    # create new transform from old transform
    new_transform = ln.Transform(name="My new transform", is_new_version_of=transform)
    assert transform.version == "1"
    assert new_transform.stem_uid == transform.stem_uid
    assert new_transform.version == "2"

    transform.delete()


def test_delete():
    # prepare the creation of a transform with its artifacts
    transform = ln.Transform(name="My transform")
    transform.save()
    run = ln.Run(transform)
    report_path = Path("report.html")
    with open(report_path, "w") as f:
        f.write("a")
    source_code_path = Path("code.py")
    with open(source_code_path, "w") as f:
        f.write("b")
    environment_path = Path("environment.txt")
    with open(environment_path, "w") as f:
        f.write("c")
    report = ln.Artifact(report_path, description=f"Report of {run.uid}")
    report.save()
    report_path.unlink()
    report_path = report.path
    source_code = ln.Artifact(
        source_code_path, description=f"Source of {transform.uid}"
    )
    source_code.save()
    source_code_path.unlink()
    source_code_path = source_code.path
    environment = ln.Artifact(environment_path, description="requirement.txt")
    environment.save()
    environment_path.unlink()
    environment_path = environment.path
    transform.latest_report = report
    transform.source_code = source_code
    transform.save()
    run.report = report
    run.environment = environment
    run.save()
    assert report_path.exists()
    assert source_code_path.exists()
    assert environment_path.exists()
    # now delete everything
    transform.delete()
    assert not report_path.exists()
    assert not source_code_path.exists()
    assert not environment_path.exists()
    assert (
        len(
            ln.Artifact.filter(id__in=[report.id, source_code.id, environment.id]).all()
        )
        == 0
    )
    assert len(ln.Run.filter(id=run.id).all()) == 0
