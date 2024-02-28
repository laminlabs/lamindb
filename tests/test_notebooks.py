from pathlib import Path

import nbproject_test

notebook_dir = Path(__file__).parent / "notebooks/"


def test_no_title():
    nbproject_test.execute_notebooks(notebook_dir / "no-title.ipynb")
