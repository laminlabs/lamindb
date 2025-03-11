import re

import lamindb as ln
import pytest
from lamindb.errors import FieldValidationError


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

    with pytest.raises(
        ValueError,
        match=re.escape(
            "'my_type' should start with a capital letter given you're defining a type"
        ),
    ):
        ln.ULabel(name="my_type", is_type=True)


def test_ulabel_plural_type_warning(ccaplog):
    ln.ULabel(name="MyThings", is_type=True)
    assert (
        "name 'MyThings' for type ends with 's', in case you're naming with plural, consider the singular for a type name"
        in ccaplog.text
    )
