import sys
from pathlib import Path

import nbproject_test as test

sys.path[:0] = [str(Path(__file__).parent.parent)]

from noxfile import GROUPS

DOCS = Path(__file__).parents[1] / "docs/"


def test_tutorial():
    for artifactname in GROUPS["tutorial"]:
        test.execute_notebooks(DOCS / artifactname, write=True)


def test_guide():
    for artifactname in GROUPS["guide"]:
        test.execute_notebooks(DOCS / artifactname, write=True)


def test_biology():
    for artifactname in GROUPS["biology"]:
        test.execute_notebooks(DOCS / artifactname, write=True)
