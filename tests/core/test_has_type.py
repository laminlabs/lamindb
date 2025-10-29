import lamindb as ln
import pytest
from django.db import IntegrityError


def test_invalid_type_record():
    # also see test_invalid_type_record_with_schema in test_record.py
    no_record_type = ln.Record(name="no_type").save()
    with pytest.raises(ValueError) as error:
        ln.Record(name="with_invalid_type", type=no_record_type).save()
    assert error.exconly().startswith(
        "ValueError: You can only assign a record of `is_type=True` as `type` to another record"
    )
    # test at the database level
    no_record_type.is_type = True
    with pytest.raises(IntegrityError) as error:
        ln.Record(name="with_invalid_type", type=no_record_type).save()
    assert "record_type_is_valid_fk" in error.exconly()
    no_record_type.delete(permanent=True)


def test_invalid_type_feature():
    no_feature_type = ln.Feature(name="no_type", dtype="str").save()
    with pytest.raises(IntegrityError) as error:
        ln.Feature(name="with_invalid_type", dtype="str", type=no_feature_type).save()
    assert "feature_type_is_valid_fk" in error.exconly()
    no_feature_type.delete(permanent=True)


def test_invalid_type_schema():
    no_schema_type = ln.Schema(name="no_type", itype=ln.Feature).save()
    with pytest.raises(IntegrityError) as error:
        ln.Schema(
            name="with_invalid_type", type=no_schema_type, itype=ln.Feature
        ).save()
    assert "schema_type_is_valid_fk" in error.exconly()
    no_schema_type.delete(permanent=True)


def test_invalid_type_project():
    no_project_type = ln.Project(name="no_type").save()
    with pytest.raises(IntegrityError) as error:
        ln.Project(name="with_invalid_type", type=no_project_type).save()
    assert "project_type_is_valid_fk" in error.exconly()
    no_project_type.delete(permanent=True)


def test_invalid_type_reference():
    no_reference_type = ln.Reference(name="no_type").save()
    with pytest.raises(IntegrityError) as error:
        ln.Reference(name="with_invalid_type", type=no_reference_type).save()
    assert "reference_type_is_valid_fk" in error.exconly()
    no_reference_type.delete(permanent=True)


def test_invalid_type_ulabel():
    no_ulabel_type = ln.ULabel(name="no_type").save()
    with pytest.raises(IntegrityError) as error:
        ln.ULabel(name="with_invalid_type", type=no_ulabel_type).save()
    assert "ulabel_type_is_valid_fk" in error.exconly()
    no_ulabel_type.delete(permanent=True)
