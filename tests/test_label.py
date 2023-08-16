import pytest

import lamindb as ln


def test_label():
    with pytest.raises(ValueError):
        ln.Label(x=1)
