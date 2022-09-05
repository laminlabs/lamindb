from pathlib import Path

from nbproject.dev import test

from lamindb._logger import logger


def test_notebooks():
    # assuming this is in the tests folder
    docs_folder = Path(__file__).parents[1] / "docs/"

    for check_folder in docs_folder.glob("./**"):
        # these are the typical notebook testpaths
        if not str(check_folder).endswith(("faq", "guide")):
            continue
        logger.info(f"\n---{check_folder.stem}---")
        test.execute_notebooks(check_folder, write=True)
