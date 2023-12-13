import lamindb as ln
import pytest


def test_ulabel():
    with pytest.raises(ValueError) as error:
        ln.ULabel(x=1)
    assert (
        error.exconly()
        == "ValueError: Only name, description, reference, reference_type are valid keyword arguments"
    )

    with pytest.raises(ValueError) as error:
        ln.ULabel(1)
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed"
