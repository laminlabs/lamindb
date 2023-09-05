import pytest

import lamindb as ln


def test_label():
    with pytest.raises(ValueError) as error:
        ln.Label(x=1)
    error.exconly() == "ValueError: Only name & description are valid keyword arguments"

    with pytest.raises(ValueError) as error:
        ln.Label(1)
    error.exconly() == "ValueError: Only one non-keyword arg allowed"
