from pathlib import Path

import nbproject_test as test
from lamin_logger import logger


def test_notebooks():
    # assuming this is in the tests folder
    docsdir = Path(__file__).parents[1] / "docs/"

    for subdir in ["guide", "faq"]:
        checkdir = docsdir / subdir
        logger.info(f"\n---{checkdir.stem}---")
        test.execute_notebooks(checkdir, write=True)
