import os
import subprocess
from pathlib import Path

import nbproject_test

notebook_dir = "./sub/lamin-cli/tests/notebooks/"


def test_track_not_initialized():
    env = os.environ
    env["LAMIN_TESTING"] = "true"
    result = subprocess.run(
        f"lamin track {notebook_dir}not-initialized.ipynb",
        shell=True,
        capture_output=True,
        env=env,
    )
    assert result.returncode == 0
    assert "attached notebook id to ipynb file" in result.stdout.decode()


def test_track_no_title():
    env = os.environ
    env["LAMIN_TESTING"] = "true"
    result = subprocess.run(
        f"lamin track {notebook_dir}no-title.ipynb",
        shell=True,
        capture_output=True,
        env=env,
    )
    assert result.returncode == 0
    # see lamindb_setup/_notebooks.py::update_notebook_metadata
    assert "updated notebook metadata" in result.stdout.decode()


def test_save_no_title():
    env = os.environ
    env["LAMIN_TESTING"] = "true"
    result = subprocess.run(
        f"lamin save {notebook_dir}no-title.ipynb",
        shell=True,
        capture_output=True,
        env=env,
    )
    assert result.returncode == 1
    assert "No title!" in result.stdout.decode()


def test_save_non_consecutive():
    env = os.environ
    env["LAMIN_TESTING"] = "true"
    result = subprocess.run(
        f"lamin save {notebook_dir}with-title-and-initialized-non-consecutive.ipynb",  # noqa
        shell=True,
        capture_output=True,
        env=env,
    )
    assert result.returncode == 1
    assert "Aborted (non-consecutive)!" in result.stdout.decode()


def test_save_consecutive():
    notebook_path = Path(
        f"{notebook_dir}with-title-and-initialized-consecutive.ipynb"
    ).resolve()
    env = os.environ
    env["LAMIN_TESTING"] = "true"

    # let's try to save a notebook for which `ln.track()` was never run
    # it's going to fail
    result = subprocess.run(
        f"lamin save {str(notebook_path)}",
        shell=True,
        capture_output=True,
        env=env,
    )
    assert result.returncode == 1
    assert "Didn't find notebook with stem_id" in result.stdout.decode()

    # now, let's re-run this notebook so that ln.track() is actually run
    nbproject_test.execute_notebooks(notebook_path)
    # and save again
    result = subprocess.run(
        f"lamin save {str(notebook_path)}",
        shell=True,
        capture_output=True,
        env=env,
    )
    print(result.stdout)
    print(result.stderr)
    assert result.returncode == 0
    assert (
        "saved notebook and wrote source file and html report" in result.stdout.decode()
    )

    # now, assume the user modifies the notebook and saves
    # it without changing id or version
    from nbproject.dev import read_notebook, write_notebook

    nb = read_notebook(notebook_path)
    # duplicate last cell
    new_cell = nb.cells[-1].copy()
    new_cell["execution_count"] += 1
    nb.cells.append(new_cell)
    write_notebook(nb, notebook_path)
    result = subprocess.run(
        f"lamin save {str(notebook_path)}",
        shell=True,
        capture_output=True,
        env=env,
    )
    print(result.stdout)
    print(result.stderr)
    assert result.returncode == 0
    assert (
        "saved notebook and wrote source file and html report" in result.stdout.decode()
    )
