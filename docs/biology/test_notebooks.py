from pathlib import Path

import nbproject_test as test

import lamindb as ln


def test_notebooks():
    nbdir = Path(__file__).parent
    ln.setup.login("testuser1")
    test.execute_notebooks(nbdir, write=True)


def test_lnschema_bionty():
    notebook = Path(__file__).parent.parent / "lnschema-bionty.ipynb"
    test.execute_notebooks(notebook, write=True)
