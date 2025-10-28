import re

import lamindb as ln
import pytest
from django.db import IntegrityError
from lamindb.errors import FieldValidationError


def test_record():
    with pytest.raises(
        FieldValidationError,
        match=re.escape(
            "Only name, type, is_type, description, schema, reference, reference_type are valid keyword arguments"
        ),
    ):
        ln.Record(x=1)

    with pytest.raises(ValueError) as error:
        ln.Record(1)
    assert error.exconly() == "ValueError: Only one non-keyword arg allowed"

    with pytest.raises(
        ValueError,
        match=re.escape(
            "'my_type' should start with a capital letter given you're defining a type"
        ),
    ):
        ln.Record(name="my_type", is_type=True)


def test_record_plural_type_warning(ccaplog):
    ln.Record(name="MyThings", is_type=True)
    assert (
        "name 'MyThings' for type ends with 's', in case you're naming with plural, consider the singular for a type name"
        in ccaplog.text
    )


def test_name_lookup():
    my_type = ln.Record(name="MyType", is_type=True).save()
    label1 = ln.Record(name="label 1", type=my_type).save()
    label2 = ln.Record(name="label 1", type=my_type)
    assert label2 == label1
    label2 = ln.Record(name="label 1")
    assert label2 != label1
    label2.save()
    label3 = ln.Record(name="label 1")
    assert label3 == label2
    label2.delete(permanent=True)
    label1.delete(permanent=True)
    my_type.delete(permanent=True)


def test_invalid_type_record_with_schema():
    schema = ln.Schema(name="test_schema", itype=ln.Feature).save()

    record_type_with_schema = ln.Record(
        name="TypeWithSchema", is_type=True, schema=schema
    ).save()

    with pytest.raises(IntegrityError) as error:
        ln.Record(name="InvalidType", is_type=True, type=record_type_with_schema).save()
    assert "record_type_is_valid_fk" in error.exconly()

    record_type_with_schema.delete(permanent=True)
    schema.delete(permanent=True)
