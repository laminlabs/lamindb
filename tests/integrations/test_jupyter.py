import builtins

import pytest
from lamindb.integrations.jupyter import read_notebook


def test_read_notebook_missing_deps(monkeypatch):
    real_import = builtins.__import__

    def fake_import(name, *args, **kwargs):
        if name in {"jupytext", "nbformat"}:
            raise ImportError(name)
        return real_import(name, *args, **kwargs)

    monkeypatch.setattr(builtins, "__import__", fake_import)
    with pytest.raises(ImportError, match="install nbconvert & jupytext"):
        read_notebook("notebook.ipynb")
