import pytest

import lamindb as ln


def test_create():
    transform = ln.Transform(stem_id="NJvdsWWbJlZS", version="0")
    assert transform.id == "NJvdsWWbJlZSz8"
    transform = ln.Transform(id="wNWf6CCdoafkz8")
    assert transform.stem_id == "wNWf6CCdoafk"
    with pytest.raises(ValueError):
        ln.Transform(version=1)
