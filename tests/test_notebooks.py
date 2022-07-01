from pathlib import Path

from nbproject.dev import test


def test_notebooks():
    # assuming this is in the tests folder
    docs_folder = Path(__file__).parents[1] / "docs/"

    for check_folder in docs_folder.glob("./*"):
        if check_folder.is_dir():
            test.execute_notebooks(check_folder, write=True)
