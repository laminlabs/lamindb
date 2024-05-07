import subprocess

import lamindb as ln
import pytest
from lamindb.core._run_context import get_uid_ext, run_context


def test_track_with_multi_parents():
    parent1 = ln.Transform(name="parent 1")
    parent1.save()
    parent2 = ln.Transform(name="parent 2")
    parent2.save()
    child = ln.Transform(name="Child")
    child.save()
    child.parents.set([parent1, parent2])
    params = {"param1": 1, "param2": "my-string", "param3": 3.14}
    ln.track(transform=child, params=params)

    assert ln.core.run_context.run.json == params
    # unset to remove side effects
    ln.core.run_context.run = None
    ln.core.run_context.transform = None


def test_track_notebook_colab():
    notebook_path = "/fileId=1KskciVXleoTeS_OGoJasXZJreDU9La_l"
    ln.core.run_context._track_notebook(path=notebook_path)


def test_create_or_load_transform(monkeypatch):
    title = "title"
    stem_uid = "NJvdsWWbJlZS"
    version = "0"
    uid = "NJvdsWWbJlZS6K79"
    assert uid == f"{stem_uid}{get_uid_ext(version)}"
    run_context._create_or_load_transform(
        stem_uid=stem_uid,
        version=version,
        name=title,
        transform_type="notebook",
    )
    assert run_context.transform.uid == uid
    assert run_context.transform.version == version
    assert run_context.transform.name == title
    run_context._create_or_load_transform(
        transform=run_context.transform,
        stem_uid=stem_uid,
        version=version,
        name=title,
    )
    assert run_context.transform.uid == uid
    assert run_context.transform.version == version
    assert run_context.transform.name == title

    # now, test an updated transform name (updated notebook title)

    # monkeypatch the "input" function, so that it returns "n"
    # this simulates the user entering "n" in the terminal
    monkeypatch.setattr("builtins.input", lambda _: "n")
    run_context._create_or_load_transform(
        transform=run_context.transform,
        stem_uid=stem_uid,
        version=version,
        name="updated title",
    )
    assert run_context.transform.uid == uid
    assert run_context.transform.version == version
    assert run_context.transform.name == "updated title"

    # test the user responding with "y"
    monkeypatch.setattr("builtins.input", lambda _: "y")
    with pytest.raises(SystemExit) as error:
        run_context._create_or_load_transform(
            transform=run_context.transform,
            stem_uid=stem_uid,
            version=version,
            name="updated title again",
        )
    assert (
        "UpdateTransformSettings: Please update your transform settings as follows"
        in error.exconly()
    )


def test_run_script():
    script_path = "sub/lamin-cli/tests/scripts/run-track-and-finish-sync-git.py"
    result = subprocess.run(
        f"python {script_path}",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 0
    assert "saved: Transform" in result.stdout.decode()
    assert "saved: Run" in result.stdout.decode()
    transform = ln.Transform.filter(name="run-track-and-finish-sync-git.py").one()
    # the algorithm currently picks different commits depending on the state of the repo
    # any of these commits are valid
    assert transform.uid == "m5uCHTTpJnjQ5zKv"
    assert transform.reference.endswith(
        "/tests/scripts/run-track-and-finish-sync-git.py"
    )
    assert transform.reference.startswith(
        "https://github.com/laminlabs/lamin-cli/blob/"
    )
    assert transform.reference_type == "url"
    # ensure that the source code is not saved as an output artifact
    assert transform.latest_run.output_artifacts.count() == 0
    assert transform.runs.count() == 1
    assert transform.source_code.hash == "-QN2dVdC8T3xWG8vBl-wew"
    assert transform.source_code.transform is None
    assert transform.source_code.run is None


@pytest.mark.parametrize("type", ["notebook", "script"])
def test_track_notebook_or_script_manually(type):
    transform = ln.Transform(name="My notebook", type=type)
    with pytest.raises(ValueError) as error:
        ln.track(transform=transform)
    assert (
        error.exconly()
        == "ValueError: Use ln.track() without passing transform in a notebook or script - metadata is automatically parsed"
    )
