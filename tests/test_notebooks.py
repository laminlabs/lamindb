from pathlib import Path

import nbproject_test as test
from lamin_logger import logger

import lamindb as ln


def test_notebooks():
    # assuming this is in the tests folder
    docsdir = Path(__file__).parents[1] / "docs/"

    ln.setup.init(storage=docsdir / "guide" / "mydata", schema="bionty,lamin1")

    for subdir in ["guide", "faq"]:
        checkdir = docsdir / subdir
        logger.info(f"\n---{checkdir.stem}---")
        test.execute_notebooks(checkdir, write=True)
