import re

import bionty as bt
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


def test_record_features_add_remove_values():
    record_type1 = ln.Record(name="RecordType1", is_type=True).save()
    record_entity1 = ln.Record(name="entity1", type=record_type1).save()
    record_entity2 = ln.Record(name="entity2", type=record_type1).save()
    artifact = ln.Artifact(".gitignore", key="test-artifact").save()
    transform = ln.Transform(name="test-transform").save()
    run = ln.Run(transform, name="test-run").save()

    feature_str = ln.Feature(name="feature_str", dtype=str).save()
    feature_int = ln.Feature(name="feature_int", dtype=int).save()
    feature_dict = ln.Feature(name="feature_dict", dtype=dict).save()
    feature_type1 = ln.Feature(name="feature_type1", dtype=record_type1).save()
    feature_type2 = ln.Feature(name="feature_type2", dtype=list[record_type1]).save()
    feature_user = ln.Feature(name="feature_user", dtype=ln.User).save()
    feature_project = ln.Feature(name="feature_project", dtype=ln.Project).save()
    feature_artifact = ln.Feature(name="feature_artifact", dtype=ln.Artifact).save()
    feature_run = ln.Feature(name="feature_run", dtype=ln.Run.uid).save()
    feature_cell_line = ln.Feature(name="feature_cell_line", dtype=bt.CellLine).save()
    feature_cl_ontology_id = ln.Feature(
        name="feature_cl_ontology_id", dtype=bt.CellLine.ontology_id
    ).save()

    test_record = ln.Record(name="test_record").save()
    test_project = ln.Project(name="test_project").save()
    hek293 = bt.CellLine.from_source(name="HEK293").save()

    # no schema validation

    test_values = {
        "feature_str": "a string value",
        "feature_int": 42,
        "feature_dict": {"key": "value", "number": 123, "list": [1, 2, 3]},
        "feature_type1": "entity1",
        "feature_type2": ["entity1", "entity2"],
        "feature_user": ln.setup.settings.user.handle,
        "feature_project": "test_project",
        "feature_cell_line": "HEK293",
        "feature_cl_ontology_id": "CLO:0001230",
        "feature_artifact": "test-artifact",
        "feature_run": run.uid,
    }

    test_record.features.add_values(test_values)
    assert test_record.features.get_values() == test_values

    # remove values

    test_record.features.remove_values("feature_int")
    test_values.pop("feature_int")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_type1")
    test_values.pop("feature_type1")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_type2")
    test_values.pop("feature_type2")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_cell_line")
    test_values.pop("feature_cell_line")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_user")
    test_values.pop("feature_user")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_artifact")
    test_values.pop("feature_artifact")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_run")
    test_values.pop("feature_run")
    assert test_record.features.get_values() == test_values

    # schema validation

    feature_str = ln.Feature.get(name="feature_str")
    feature_int = ln.Feature.get(name="feature_int")
    schema = ln.Schema([feature_str, feature_int], name="test_schema").save()
    test_form = ln.Record(name="TestForm", is_type=True, schema=schema).save()
    test_record_in_form = ln.Record(name="test_record_in_form", type=test_form).save()
    with pytest.raises(ln.errors.ValidationError) as error:
        test_record_in_form.features.add_values({"feature_type1": "entity1"})
    assert "COLUMN_NOT_IN_DATAFRAME" in error.exconly()

    test_record_in_form.delete(permanent=True)
    test_form.delete(permanent=True)
    schema.delete(permanent=True)
    test_record.delete(permanent=True)
    feature_str.delete(permanent=True)
    feature_int.delete(permanent=True)
    feature_type1.delete(permanent=True)
    feature_type2.delete(permanent=True)
    feature_user.delete(permanent=True)
    feature_project.delete(permanent=True)
    feature_dict.delete(permanent=True)
    feature_artifact.delete(permanent=True)
    feature_run.delete(permanent=True)
    record_entity1.delete(permanent=True)
    record_entity2.delete(permanent=True)
    record_type1.delete(permanent=True)
    test_project.delete(permanent=True)
    feature_cell_line.delete(permanent=True)
    feature_cl_ontology_id.delete(permanent=True)
    hek293.delete(permanent=True)
    artifact.delete(permanent=True)
    run.delete(permanent=True)
    transform.delete(permanent=True)
