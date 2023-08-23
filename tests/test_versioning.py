import pytest

from lamindb.dev.versioning import set_version


def test_set_version():
    # all remaining lines are covered in notebooks
    with pytest.raises(ValueError):
        set_version(None, "1.2")
    assert set_version(None, "0") == "1"
    assert set_version(None, "1") == "2"
    assert set_version("1.2.3", "0") == "1.2.3"
    assert set_version("1.2.3") == "1.2.3"
