import os
from pathlib import Path

import nbproject_test

notebook_dir = Path(__file__).parent / "notebooks/"


def test_all_notebooks():
    env = os.environ
    env["LAMIN_TESTING"] = "true"
    nbproject_test.execute_notebooks(notebook_dir)
    del env["LAMIN_TESTING"]
