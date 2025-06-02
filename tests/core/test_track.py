import subprocess
import sys
import time
from pathlib import Path

import lamindb as ln
import lamindb_setup as ln_setup
import pytest
from lamindb._finish import clean_r_notebook_html, get_shortcut
from lamindb.core._context import LogStreamTracker, context
from lamindb.errors import TrackNotCalled, ValidationError

SCRIPTS_DIR = Path(__file__).parent.resolve() / "scripts"
NOTEBOOKS_DIR = Path(__file__).parent.resolve() / "notebooks"


def test_track_basic_invocation():
    project = "non-existing project"
    with pytest.raises(ln.errors.InvalidArgument) as error:
        ln.track(project=project)
    assert (
        error.exconly()
        == f"lamindb.errors.InvalidArgument: Project '{project}' not found, either create it with `ln.Project(name='...').save()` or fix typos."
    )
    space = "non-existing space"
    with pytest.raises(ln.errors.InvalidArgument) as error:
        ln.track(space=space)
    assert (
        error.exconly()
        == f"lamindb.errors.InvalidArgument: Space '{space}', please check on the hub UI whether you have the correct `uid` or `name`."
    )

    predecessor1 = ln.Transform(key="parent 1").save()
    predecessor2 = ln.Transform(key="parent 2").save()
    successor = ln.Transform(key="successor").save()
    successor.predecessors.set([predecessor1, predecessor2])

    # first invocation
    params = {"param1": 1, "param2": "my-string", "param3": 3.14}
    with pytest.raises(ValidationError) as exc:
        ln.track(transform=successor, params=params)
    assert (
        exc.exconly()
        == """lamindb.errors.ValidationError: These keys could not be validated: ['param1', 'param2', 'param3']
Here is how to create a feature:

  ln.Feature(name='param1', dtype='int').save()
  ln.Feature(name='param2', dtype='cat ? str').save()
  ln.Feature(name='param3', dtype='float').save()"""
    )
    ln.Feature(name="param1", dtype="int").save()
    ln.Feature(name="param2", dtype="str").save()
    ln.Feature(name="param3", dtype="float").save()
    ln.track(transform=successor, params=params)
    print("outside", id(ln.context))
    assert ln.context.run.features.get_values() == params
    # second invocation
    params = {"param1": 1, "param2": "my-string", "param3": 3.14, "param4": [1, 2]}
    param4 = ln.Feature(name="param4", dtype="int").save()
    with pytest.raises(ValidationError) as exc:
        ln.track(transform=successor, params=params)
    assert (
        exc.exconly()
        == """lamindb.errors.ValidationError: Expected dtype for 'param4' is 'int', got 'list[int]'"""
    )
    # fix param4 dtype
    param4.dtype = "list[int]"
    param4.save()
    # re-run
    ln.track(transform=successor, params=params)
    assert ln.context.run.features.get_values() == params

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
    ln.context._track_notebook(path_str=notebook_path)


def test_finish_before_track():
    ln.context._run = None
    with pytest.raises(TrackNotCalled) as error:
        ln.finish()
    assert "Please run `ln.track()` before `ln.finish()" in error.exconly()


def test_invalid_transform_type():
    transform = ln.Transform(key="test transform")
    ln.track(transform=transform)
    ln.context._path = None
    ln.context.run.transform.type = "script"
    with pytest.raises(ValueError) as error:
        ln.finish()
    assert "Transform type is not allowed to be" in error.exconly()

    # unset to remove side effects
    ln.context._run = None


def test_create_or_load_transform():
    title = "title"
    version = "2.0"
    uid = "NJvdsWWbJlZS0000"
    context.uid = uid
    context.version = version
    context._path = Path("my-test-transform-create-or-load.py")
    context._path.touch(exist_ok=True)
    context._create_or_load_transform(
        description=title,
        transform_type="notebook",
    )
    assert context._transform.uid == uid
    assert context._transform.version == version
    assert context._transform.description == title
    context._create_or_load_transform(
        description=title,
    )
    assert context._transform.uid == uid
    assert context._transform.version == version
    assert context._transform.description == title

    # now, test an updated transform name
    context._create_or_load_transform(
        description="updated title",
    )
    assert context._transform.uid == uid
    assert context._transform.version == version
    assert context._transform.description == "updated title"

    # unset to remove side effects
    ln.context._uid = None
    ln.context._run = None
    ln.context._transform = None
    ln.context._path.unlink()
    ln.context._path = None


def test_run_scripts():
    # regular execution
    result = subprocess.run(  # noqa: S602
        f"python {SCRIPTS_DIR / 'script-to-test-versioning.py'}",
        shell=True,
        capture_output=True,
    )
    assert result.returncode == 0
    assert (
        "created Transform('Ro1gl7n8YrdH0000'), started new Run("
        in result.stdout.decode()
    )

    # updated key (filename change)
    result = subprocess.run(  # noqa: S602
        f"python {SCRIPTS_DIR / 'script-to-test-filename-change.py'}",
        shell=True,
        capture_output=True,
    )
    assert result.returncode == 0
    assert "renaming transform" in result.stdout.decode()

    # version already taken
    result = subprocess.run(  # noqa: S602
        f"python {SCRIPTS_DIR / 'duplicate1/script-to-test-versioning.py'}",
        shell=True,
        capture_output=True,
    )
    assert result.returncode == 1
    assert (
        "✗ version '1' is already taken by Transform('Ro1gl7n8YrdH0000'); please set another version, e.g., ln.context.version = '1.1'"
        in result.stderr.decode()
    )

    # regular version bump
    result = subprocess.run(  # noqa: S602
        f"python {SCRIPTS_DIR / 'duplicate2/script-to-test-versioning.py'}",
        shell=True,
        capture_output=True,
    )
    assert result.returncode == 0
    assert (
        "created Transform('Ro1gl7n8YrdH0001'), started new Run("
        in result.stdout.decode()
    )
    assert not ln.Transform.get("Ro1gl7n8YrdH0000").is_latest
    assert ln.Transform.get("Ro1gl7n8YrdH0001").is_latest

    # inconsistent version
    result = subprocess.run(  # noqa: S602
        f"python {SCRIPTS_DIR / 'duplicate3/script-to-test-versioning.py'}",
        shell=True,
        capture_output=True,
    )
    assert result.returncode == 1
    assert (
        "✗ please pass consistent version: ln.context.version = '2'"
        in result.stderr.decode()
    )

    # multiple folders, no match
    ln.Transform.filter(key__endswith="script-to-test-versioning.py").update(
        key="teamA/script-to-test-versioning.py"
    )
    # this test creates a transform with key script-to-test-versioning.py at the root level
    result = subprocess.run(  # noqa: S602
        f"python {SCRIPTS_DIR / 'duplicate4/script-to-test-versioning.py'}",
        shell=True,
        capture_output=True,
    )
    assert result.returncode == 0
    assert "ignoring transform" in result.stdout.decode()

    # multiple folders, match
    transform = ln.Transform(
        description="Dummy title",
        key="duplicate4/script-to-test-versioning.py",
        type="script",
    ).save()
    result = subprocess.run(  # noqa: S602
        f"python {SCRIPTS_DIR / 'duplicate4/script-to-test-versioning.py'}",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 0
    assert f"{transform.stem_uid}" in result.stdout.decode()
    assert "making new version" not in result.stdout.decode()

    transform.source_code = "dummy"
    transform.save()

    # multiple folders, match and transform has saved source code
    result = subprocess.run(  # noqa: S602
        f"python {SCRIPTS_DIR / 'duplicate4/script-to-test-versioning.py'}",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 0
    assert f"{transform.stem_uid}" in result.stdout.decode()
    assert "making new version" in result.stdout.decode()


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
    assert "started new Run" in result.stdout.decode()
    transform = ln.Transform.get(key="run-track-and-finish-sync-git.py")
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
    assert transform.description == "My good script"
    # ensure that the source code is not saved as an output artifact
    assert transform.latest_run.output_artifacts.count() == 0
    assert transform.runs.count() == 1
    assert transform.hash == "VC1oTPcaVSrzNrXUT9p4qw"


@pytest.mark.parametrize("type", ["notebook", "script"])
def test_track_notebook_or_script_manually(type):
    transform = ln.Transform(key="My notebook", type=type)
    with pytest.raises(ValueError) as error:
        ln.track(transform=transform)
    assert (
        error.exconly()
        == "ValueError: Use `ln.track()` without passing transform in a notebook or script - metadata is automatically parsed"
    )


def test_clean_r_notebook_html():
    orig_notebook_path = NOTEBOOKS_DIR / "basic-r-notebook.Rmd.html"
    content = orig_notebook_path.read_text()
    orig_notebook_path.write_text(content.replace("SHORTCUT", get_shortcut()))
    comparison_path = NOTEBOOKS_DIR / "basic-r-notebook.Rmd.cleaned.html"
    compare = comparison_path.read_text()
    comparison_path.unlink()
    title_text, cleaned_path = clean_r_notebook_html(orig_notebook_path)
    assert comparison_path == cleaned_path
    assert title_text == "My exemplary R analysis"
    assert compare == cleaned_path.read_text()  # check that things have been stripped
    comparison_path.write_text(compare)
    orig_notebook_path.write_text(content.replace(get_shortcut(), "SHORTCUT"))


class MockRun:
    def __init__(self, uid):
        self.uid = uid
        self.report = None
        self.saved = False

    def save(self):
        self.saved = True


def test_logstream_tracker_multiple():
    tracker1 = LogStreamTracker()
    tracker2 = LogStreamTracker()
    tracker3 = LogStreamTracker()

    try:
        # Start trackers one by one and print messages
        print("Initial stdout")

        tracker1.start(MockRun("run1"))
        print("After starting tracker1")

        tracker2.start(MockRun("run2"))
        print("After starting tracker2")

        tracker3.start(MockRun("run3"))
        print("After starting tracker3")

        print("Testing stderr", file=sys.stderr)

        time.sleep(0.1)

        # Clean up in reverse order
        tracker3.finish()
        tracker2.finish()
        tracker1.finish()

        # Verify log contents - each log should only contain messages after its start
        expected_contents = {
            1: [
                "After starting tracker1",
                "After starting tracker2",
                "After starting tracker3",
                "Testing stderr",
            ],
            2: ["After starting tracker2", "After starting tracker3", "Testing stderr"],
            3: ["After starting tracker3", "Testing stderr"],
        }

        for i in range(1, 4):
            log_path = Path(ln_setup.settings.cache_dir / f"run_logs_run{i}.txt")
            with open(log_path) as f:
                content = f.read()
                print(f"\nContents of run{i} log:")
                print(content)
                # Check each expected line is in the content
                for expected_line in expected_contents[i]:
                    assert expected_line in content, (
                        f"Expected '{expected_line}' in log {i}"
                    )

                # Check earlier messages are NOT in the content
                if i > 1:
                    assert "Initial stdout" not in content
                    assert "After starting tracker" + str(i - 1) not in content

    finally:
        # Cleanup
        for i in range(1, 4):
            log_path = Path(ln_setup.settings.cache_dir / f"run_logs_run{i}.txt")
            if log_path.exists():
                log_path.unlink()


def test_logstream_tracker_exception_handling():
    tracker = LogStreamTracker()
    original_excepthook = sys.excepthook
    run = MockRun("error")

    try:
        tracker.start(run)
        print("Before error")

        # Create and capture exception info
        exc_type = ValueError
        exc_value = ValueError("Test error")
        exc_traceback = None
        try:
            raise exc_value
        except ValueError:
            exc_traceback = sys.exc_info()[2]

        # Handle the exception - this will trigger cleanup
        tracker.handle_exception(exc_type, exc_value, exc_traceback)

        # Verify run status
        assert run.saved
        assert run.report is not None

        # Verify the content was written before cleanup
        content = run.report.cache().read_text()
        print("Log contents:", content)
        assert "Before error" in content
        assert "ValueError: Test error" in content
        assert "Traceback" in content

    finally:
        sys.excepthook = original_excepthook
