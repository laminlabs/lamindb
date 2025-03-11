import os
from pathlib import Path

import nbproject_test

notebook_dir = Path(__file__).parent / "notebooks/"
notebook_dir_duplicate = Path(__file__).parent / "notebooks/duplicate/"


def test_all_notebooks():
    env = os.environ
    env["LAMIN_TESTING"] = "true"
    nbproject_test.execute_notebooks(notebook_dir)
    # re-run one of the notebooks to trigger hash lookup
    nbproject_test.execute_notebooks(
        notebook_dir / "with-title-initialized-consecutive-finish-not-last-cell.ipynb"
    )
    nbproject_test.execute_notebooks(notebook_dir_duplicate)
    del env["LAMIN_TESTING"]
