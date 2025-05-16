import os
import subprocess
from pathlib import Path

import lamindb as ln
import nbproject_test

notebook_dir = Path(__file__).parent / "notebooks/"
notebook_dir_duplicate = Path(__file__).parent / "notebooks/duplicate/"


def test_all_notebooks():
    env = os.environ
    env["LAMIN_TESTING"] = "true"
    nbproject_test.execute_notebooks(notebook_dir)
    nbproject_test.execute_notebooks(notebook_dir_duplicate)
    del env["LAMIN_TESTING"]


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
        key="with-title-initialized-consecutive-finish-not-last-cell.ipynb"
    ).uid

    # now, assume the user renames the notebook
    new_path = notebook_path.with_name("no-uid-renamed.ipynb")
    os.system(f"cp {notebook_path} {new_path}")  # noqa: S605

    env = os.environ
    env["LAMIN_TESTING"] = "true"
    result = subprocess.run(  # noqa: S602
        f"jupyter nbconvert --to notebook --inplace --execute {new_path}",
        shell=True,
        capture_output=True,
        env=env,
    )
    print(result.stdout.decode())
    print(result.stderr.decode())
    assert result.returncode == 0
    del env["LAMIN_TESTING"]

    assert ln.Transform.get(key="no-uid-renamed.ipynb").uid == uid

    # new_path.unlink()
