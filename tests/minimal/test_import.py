import importlib.util

import pytest


def test_import_lamindb() -> None:
    import lamindb  # noqa: F401


def test_lance_extra_missing() -> None:
    if importlib.util.find_spec("lancedb") is not None:
        pytest.skip("lancedb is installed")

    from lamindb.core.storage._lance import import_lancedb

    with pytest.raises(ImportError, match=r"lamindb\[lance\]"):
        import_lancedb()
