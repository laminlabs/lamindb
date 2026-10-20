import os
import subprocess
from pathlib import Path

import lamindb as ln
import nbproject_test
import pytest

notebook_dir = Path(__file__).parent / "notebooks/"
notebook_dir_duplicate = Path(__file__).parent / "notebooks/duplicate/"


def test_all_notebooks():
    nbproject_test.execute_notebooks(notebook_dir)
    nbproject_test.execute_notebooks(notebook_dir_duplicate)


def test_run_after_rename_no_uid():
    notebook_path = (
        notebook_dir / "with-title-initialized-consecutive-finish-not-last-cell.ipynb"
    )
    result = subprocess.run(  # noqa: S602
        f"jupyter nbconvert --to notebook --inplace --execute {notebook_path}",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 0

    uid = ln.Transform.get(
        key__endswith="with-title-initialized-consecutive-finish-not-last-cell.ipynb"
    ).uid

    # now, assume the user renames the notebook
    new_path = notebook_path.with_name("no-uid-renamed.ipynb")
    os.system(f"cp {notebook_path} {new_path}")  # noqa: S605

    result = subprocess.run(  # noqa: S602
        f"jupyter nbconvert --to notebook --inplace --execute {new_path}",
        shell=True,
        capture_output=True,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 0

    assert ln.Transform.get(key__endswith="no-uid-renamed.ipynb").uid == uid

    # new_path.unlink()


@pytest.fixture
def test_run_in_marimo(monkeypatch):
    """Pretend we're inside a marimo notebook with filename='/tmp/test.py'."""
    from types import SimpleNamespace

    monkeypatch.setattr(
        "lamindb.core._context._is_running_in_marimo",
        lambda: True,
    )
    monkeypatch.setattr(
        "marimo._runtime.context.safe_get_context",
        lambda: SimpleNamespace(
            filename="test.py",
            marimo_config={"display": {"dataframes": "plain"}},
        ),
    )


def test_helper_returns_filename(test_run_in_marimo, monkeypatch, ccaplog):
    from types import SimpleNamespace

    import marimo._runtime as marimo_runtime
    from lamindb.core._context import get_marimo_notebook_path, get_notebook_path

    assert get_marimo_notebook_path() == "test.py"
    assert get_notebook_path() == (Path("test.py"), "marimo")

    monkeypatch.setattr(
        "marimo._runtime.context.safe_get_context",
        lambda: SimpleNamespace(
            filename="rich.py",
            marimo_config={"display": {"dataframes": "rich"}},
        ),
    )
    assert get_marimo_notebook_path() == "rich.py"
    assert "Dataframe viewer" in ccaplog.text

    monkeypatch.setattr("marimo._runtime.context.safe_get_context", lambda: None)
    assert get_marimo_notebook_path() is None

    monkeypatch.setattr(marimo_runtime, "context", SimpleNamespace())
    with pytest.raises(RuntimeError, match="safe_get_context"):
        get_marimo_notebook_path()


def test_get_notebook_path_reports_inference_failure(monkeypatch):
    from lamindb.core._context import get_notebook_path

    def raise_error(return_env=True):
        raise OSError("no notebook")

    monkeypatch.setattr("lamindb.integrations.jupyter.notebook_path", raise_error)
    with pytest.raises(RuntimeError, match="failed to infer notebook path"):
        get_notebook_path()

    monkeypatch.setattr(
        "lamindb.integrations.jupyter.notebook_path",
        lambda return_env=True: (None, "lab"),
    )
    with pytest.raises(RuntimeError, match="failed to infer notebook path"):
        get_notebook_path()
