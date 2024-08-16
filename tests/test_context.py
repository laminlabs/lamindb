import subprocess

import lamindb as ln
import pytest
from lamindb.core._context import context, get_uid_ext
from lamindb.core.exceptions import TrackNotCalled, ValidationError


def test_track_with_multi_parents():
    parent1 = ln.Transform(name="parent 1")
    parent1.save()
    parent2 = ln.Transform(name="parent 2")
    parent2.save()
    child = ln.Transform(name="Child")
    child.save()
    child.predecessors.set([parent1, parent2])

    # first invocation
    params = {"param1": 1, "param2": "my-string", "param3": 3.14}
    with pytest.raises(ValidationError) as error:
        ln.context.track(transform=child, params=params)
    assert (
        error.exconly()
        == """lamindb.core.exceptions.ValidationError: These keys could not be validated: ['param1', 'param2', 'param3']
Here is how to create a param:

  ln.Param(name='param1', dtype='int').save()
  ln.Param(name='param2', dtype='str').save()
  ln.Param(name='param3', dtype='float').save()"""
    )
    ln.Param(name="param1", dtype="int").save()
    ln.Param(name="param2", dtype="str").save()
    ln.Param(name="param3", dtype="float").save()
    ln.context.track(transform=child, params=params)
    print("outside", id(ln.context))
    assert ln.context.run.params.get_values() == params
    # second invocation
    params = {"param1": 1, "param2": "my-string", "param3": 3.14, "param4": [1, 2]}
    param4 = ln.Param(name="param4", dtype="int").save()
    with pytest.raises(ValidationError) as error:
        ln.context.track(transform=child, params=params)
    assert (
        error.exconly()
        == """lamindb.core.exceptions.ValidationError: Expected dtype for 'param4' is 'int', got 'list[int]'"""
    )
    # fix param4 dtype
    param4.dtype = "list[int]"
    param4.save()
    # re-run
    ln.context.track(transform=child, params=params)
    assert ln.context.run.params.get_values() == params

    # test that run populates things like ULabels etc.
    ulabel = ln.ULabel(name="my-label-in-track")
    assert ulabel.run == ln.context.run

    # test that we can call ln.finish() also for pipeline-like transforms
    assert ln.context.run.finished_at is None
    ln.finish()
    assert ln.context.run.finished_at is not None

    # unset to remove side effects
    ln.context._run = None


def test_track_notebook_colab():
    notebook_path = "/fileId=1KskciVXleoTeS_OGoJasXZJreDU9La_l"
    ln.context._track_notebook(path=notebook_path)


def test_finish_before_track():
    ln.context._run = None
    with pytest.raises(TrackNotCalled) as error:
        ln.finish()
    assert "Please run `ln.context.track()` before `ln.finish()" in error.exconly()


def test_invalid_transform_type():
    transform = ln.Transform(name="test transform")
    ln.context.track(transform=transform)
    ln.context._path = None
    ln.context.run.transform.type = "script"
    with pytest.raises(ValueError) as error:
        ln.finish()
    assert "Transform type is not allowed to be" in error.exconly()

    # unset to remove side effects
    ln.context._run = None


def test_create_or_load_transform(monkeypatch):
    title = "title"
    stem_uid = "NJvdsWWbJlZS"
    version = "0"
    uid = "NJvdsWWbJlZS6K79"
    assert uid == f"{stem_uid}{get_uid_ext(version)}"
    context._create_or_load_transform(
        uid=None,
        stem_uid=stem_uid,
        version=version,
        name=title,
        transform_type="notebook",
    )
    assert context._transform.uid == uid
    assert context._transform.version == version
    assert context._transform.name == title
    context._create_or_load_transform(
        uid=None,
        transform=context._transform,
        stem_uid=stem_uid,
        version=version,
        name=title,
    )
    assert context._transform.uid == uid
    assert context._transform.version == version
    assert context._transform.name == title

    # now, test an updated transform name
    context._create_or_load_transform(
        uid=None,
        transform=context._transform,
        stem_uid=stem_uid,
        version=version,
        name="updated title",
    )
    assert context._transform.uid == uid
    assert context._transform.version == version
    assert context._transform.name == "updated title"


def test_run_scripts_for_versioning():
    # regular execution
    result = subprocess.run(  # noqa: S602
        "python ./tests/scripts/script-to-test-versioning.py",
        shell=True,
        capture_output=True,
    )
    # print(result.stdout.decode())
    assert result.returncode == 0
    assert (
        "created Transform('Ro1gl7n8YrdH0000') & created Run('"
        in result.stdout.decode()
    )

    # updated key (filename change)
    result = subprocess.run(  # noqa: S602
        "python ./tests/scripts/script-to-test-filename-change.py",
        shell=True,
        capture_output=True,
    )
    # print(result.stderr.decode())
    assert result.returncode == 1
    assert "Script filename changed." in result.stderr.decode()

    # version already taken
    result = subprocess.run(  # noqa: S602
        "python ./tests/scripts/duplicate1/script-to-test-versioning.py",
        shell=True,
        capture_output=True,
    )
    # print(result.stderr.decode())
    assert result.returncode == 1
    assert (
        "Version '1' is already taken by Transform('Ro1gl7n8YrdH0000'); please set another version, e.g., ln.context.version = '1.1'"
        in result.stderr.decode()
    )

    # regular version bump
    result = subprocess.run(  # noqa: S602
        "python ./tests/scripts/duplicate2/script-to-test-versioning.py",
        shell=True,
        capture_output=True,
    )
    # print(result.stdout.decode())
    assert result.returncode == 0
    assert (
        "created Transform('Ro1gl7n8YrdH0001') & created Run('"
        in result.stdout.decode()
    )

    # inconsistent version
    result = subprocess.run(  # noqa: S602
        "python ./tests/scripts/duplicate3/script-to-test-versioning.py",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 1
    assert (
        "Please pass consistent version: ln.context.version = '2'"
        in result.stderr.decode()
    )


def test_run_external_script():
    script_path = "sub/lamin-cli/tests/scripts/run-track-and-finish-sync-git.py"
    result = subprocess.run(  # noqa: S602
        f"python {script_path}",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 0
    assert "created Transform" in result.stdout.decode()
    assert "created Run" in result.stdout.decode()
    transform = ln.Transform.filter(key="run-track-and-finish-sync-git.py").one()
    # the algorithm currently picks different commits depending on the state of the repo
    # any of these commits are valid
    assert transform.uid == "m5uCHTTpJnjQ0000"
    assert transform.reference.endswith(
        "/tests/scripts/run-track-and-finish-sync-git.py"
    )
    assert transform.reference.startswith(
        "https://github.com/laminlabs/lamin-cli/blob/"
    )
    assert transform.reference_type == "url"
    assert transform.name == "My good script"
    # ensure that the source code is not saved as an output artifact
    assert transform.latest_run.output_artifacts.count() == 0
    assert transform.runs.count() == 1
    assert transform._source_code_artifact.hash == "Cwk0OPOyUH5nzTiU2ISlDQ"
    assert transform._source_code_artifact.transform is None
    assert transform._source_code_artifact.run is None


@pytest.mark.parametrize("type", ["notebook", "script"])
def test_track_notebook_or_script_manually(type):
    transform = ln.Transform(name="My notebook", type=type)
    with pytest.raises(ValueError) as error:
        ln.context.track(transform=transform)
    assert (
        error.exconly()
        == "ValueError: Use ln.context.track() without passing transform in a notebook or script - metadata is automatically parsed"
    )
