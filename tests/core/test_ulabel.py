import re

import lamindb as ln
import pytest
from lamindb.core.exceptions import FieldValidationError


def test_ulabel():
    with pytest.raises(
        FieldValidationError,
        match=re.escape(
            "Only name, type, is_type, description, reference, reference_type are valid keyword arguments"
        ),
    ):
        ln.ULabel(x=1)

    with pytest.raises(ValueError) as error:
        ln.ULabel(1)
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed"
