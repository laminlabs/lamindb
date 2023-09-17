import sys
from pathlib import Path

import nbproject_test as test

sys.path[:0] = [str(Path(__file__).parent.parent)]

from noxfile import GROUPS  # noqa

DOCS = Path(__file__).parents[1] / "docs/"


def test_tutorial():
    for filename in GROUPS["tutorial"]:
        test.execute_notebooks(DOCS / filename, write=True)


def test_guide():
    for filename in GROUPS["guide"]:
        test.execute_notebooks(DOCS / filename, write=True)


def test_biology():
    for filename in GROUPS["biology"]:
        test.execute_notebooks(DOCS / filename, write=True)
