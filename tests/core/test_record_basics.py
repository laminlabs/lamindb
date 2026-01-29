import os
import re
from datetime import date, datetime

import bionty as bt
import lamindb as ln
import pandas as pd
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


@pytest.mark.skipif(
    os.getenv("LAMINDB_TEST_DB_VENDOR") == "sqlite", reason="Postgres-only"
)
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


# see test_artifact_features_add_remove_query in test_artifact_external_features_annotations.py for similar test for Artifacts (populate and query by features)
def test_record_features_add_remove_values():
    record_type1 = ln.Record(name="RecordType1", is_type=True).save()
    record_entity1 = ln.Record(name="entity1", type=record_type1).save()
    record_entity2 = ln.Record(name="entity2", type=record_type1).save()
    ulabel = ln.ULabel(name="test-ulabel").save()
    artifact = ln.Artifact(".gitignore", key="test-artifact").save()
    collection = ln.Collection(artifact, key="test-collection").save()
    transform = ln.Transform(key="test-transform").save()
    run = ln.Run(transform, name="test-run").save()

    feature_bool = ln.Feature(name="feature_bool", dtype=bool).save()
    feature_str = ln.Feature(name="feature_str", dtype=str).save()
    feature_list_str = ln.Feature(name="feature_list_str", dtype=list[str]).save()
    feature_int = ln.Feature(name="feature_int", dtype=int).save()
    feature_list_int = ln.Feature(name="feature_list_int", dtype=list[int]).save()
    feature_float = ln.Feature(name="feature_float", dtype=float).save()
    feature_list_float = ln.Feature(name="feature_list_float", dtype=list[float]).save()
    feature_num = ln.Feature(name="feature_num", dtype="num").save()
    feature_list_num = ln.Feature(name="feature_list_num", dtype="list[num]").save()
    feature_datetime = ln.Feature(name="feature_datetime", dtype=datetime).save()
    feature_date = ln.Feature(name="feature_date", dtype=datetime.date).save()
    feature_dict = ln.Feature(name="feature_dict", dtype=dict).save()
    feature_type1 = ln.Feature(name="feature_type1", dtype=record_type1).save()
    feature_type1s = ln.Feature(name="feature_type1s", dtype=list[record_type1]).save()
    feature_user = ln.Feature(name="feature_user", dtype=ln.User).save()
    feature_ulabel = ln.Feature(name="feature_ulabel", dtype=ln.ULabel).save()
    feature_project = ln.Feature(name="feature_project", dtype=ln.Project).save()
    feature_artifact = ln.Feature(name="feature_artifact", dtype=ln.Artifact).save()
    feature_collection = ln.Feature(
        name="feature_collection", dtype=ln.Collection
    ).save()
    feature_run = ln.Feature(name="feature_run", dtype=ln.Run.uid).save()
    feature_cell_line = ln.Feature(name="feature_cell_line", dtype=bt.CellLine).save()
    feature_cell_lines = ln.Feature(
        name="feature_cell_lines", dtype=list[bt.CellLine]
    ).save()
    feature_cl_ontology_id = ln.Feature(
        name="feature_cl_ontology_id", dtype=bt.CellLine.ontology_id
    ).save()
    feature_gene = ln.Feature(name="feature_gene", dtype=bt.Gene).save()

    test_record = ln.Record(name="test_record").save()
    test_project = ln.Project(name="test_project").save()
    hek293 = bt.CellLine.from_source(name="HEK293").save()
    a549 = bt.CellLine.from_source(name="A-549").save()
    tmem276 = bt.Gene.from_source(symbol="Tmem276", organism="mouse").save()

    # test feature.dtype_as_object
    assert feature_bool.dtype_as_object is bool
    assert feature_str.dtype_as_object is str
    assert feature_list_str.dtype_as_object == list[str]
    assert feature_int.dtype_as_object is int
    assert feature_list_int.dtype_as_object == list[int]
    assert feature_float.dtype_as_object is float
    assert feature_list_float.dtype_as_object == list[float]
    assert feature_num.dtype_as_object is float
    assert feature_list_num.dtype_as_object == list[float]
    assert feature_datetime.dtype_as_object == datetime
    assert feature_date.dtype_as_object == date
    assert feature_dict.dtype_as_object is dict
    assert feature_type1.dtype_as_object == record_type1
    assert feature_type1s.dtype_as_object == list[record_type1]
    assert feature_user.dtype_as_object == ln.User.handle
    assert feature_ulabel.dtype_as_object == ln.ULabel.name
    assert feature_project.dtype_as_object == ln.Project.name
    assert feature_artifact.dtype_as_object == ln.Artifact.key
    assert feature_collection.dtype_as_object == ln.Collection.key
    assert feature_run.dtype_as_object == ln.Run.uid
    assert feature_cell_line.dtype_as_object == bt.CellLine.name
    assert feature_cell_lines.dtype_as_object == list[bt.CellLine.name]
    assert feature_cl_ontology_id.dtype_as_object == bt.CellLine.ontology_id
    assert feature_gene.dtype_as_object == bt.Gene.symbol

    # no schema validation
    test_values = {
        "feature_bool": True,
        "feature_str": "00810702-0006",  # this string value could be cast to datetime! don't change!
        "feature_list_str": ["a", "list", "of", "strings"],
        "feature_int": 42,
        "feature_list_int": [1, 2, 3],
        "feature_num": 3.14,
        "feature_list_num": [2.71, 3.14, 1.61],
        "feature_float": 3.14,
        "feature_list_float": [2.71, 3.14, 1.61],
        "feature_datetime": datetime(2024, 1, 1, 12, 0, 0),
        "feature_date": date(2024, 1, 1),
        "feature_dict": {"key": "value", "number": 123, "list": [1, 2, 3]},
        "feature_type1": "entity1",
        "feature_type1s": ["entity1", "entity2"],
        "feature_ulabel": "test-ulabel",
        "feature_user": ln.setup.settings.user.handle,
        "feature_project": "test_project",
        "feature_cell_line": "HEK293",
        "feature_cell_lines": ["HEK293", "A-549"],
        "feature_gene": "Tmem276",
        "feature_cl_ontology_id": "CVCL_0045",
        "feature_artifact": "test-artifact",
        "feature_collection": "test-collection",
        "feature_run": run.uid,
    }

    test_record.features.add_values(test_values)
    assert test_record.features.get_values() == test_values

    # --- Query by features (same data as above) ---
    # Equality
    assert ln.Record.filter(feature_str=test_values["feature_str"]).one() == test_record
    assert ln.Record.filter(feature_int=42).one() == test_record
    assert ln.Record.filter(feature_type1="entity1").one() == test_record
    assert ln.Record.filter(feature_cell_line="HEK293").one() == test_record
    assert (
        ln.Record.filter(feature_str=test_values["feature_str"], feature_int=42).one()
        == test_record
    )
    # Datetime and date (filter uses ISO strings as stored in JSON)
    assert ln.Record.filter(feature_datetime="2024-01-01T12:00:00").one() == test_record
    assert ln.Record.filter(feature_date="2024-01-01").one() == test_record
    # __contains (categorical)
    assert ln.Record.filter(feature_cell_line__contains="HEK").one() == test_record
    assert ln.Record.filter(feature_type1__contains="entity").one() == test_record
    # Invalid field
    with pytest.raises(ln.errors.InvalidArgument) as error:
        ln.Record.filter(feature_str_typo="x", feature_int=42).one()
    assert error.exconly().startswith(
        "lamindb.errors.InvalidArgument: You can query either by available fields:"
    )
    # DoesNotExist (no Record named "nonexistent_entity" exists)
    with pytest.raises(ln.errors.ObjectDoesNotExist) as error:
        ln.Record.filter(feature_type1="nonexistent_entity").one()
    assert "Did not find" in error.exconly()

    # Combined filter (3 keys)
    assert (
        ln.Record.filter(
            feature_str=test_values["feature_str"],
            feature_int=42,
            feature_type1="entity1",
        ).one()
        == test_record
    )
    # Bionty: filter by record
    assert ln.Record.filter(feature_cell_line=hek293).one() == test_record
    # Bionty: filter by ontology_id string
    assert ln.Record.filter(feature_cl_ontology_id="CVCL_0045").one() == test_record
    # Bionty __contains (ontology_id)
    assert (
        ln.Record.filter(feature_cl_ontology_id__contains="0045").one() == test_record
    )
    # DoesNotExist (Record not found: feature_project)
    with pytest.raises(ln.errors.ObjectDoesNotExist) as error:
        ln.Record.filter(feature_project="nonexistent_project").one()
    assert "Did not find" in error.exconly()
    # __contains returns multiple (add second record, assert, then remove)
    value_record = ln.Record(name="query_test_value_record").save()
    value_record.features.add_values({"feature_type1": "entity2"})
    assert len(ln.Record.filter(feature_type1__contains="entity")) == 2
    value_record.features.remove_values("feature_type1")
    value_record.delete(permanent=True)
    # Numeric comparators __lt, __gt (int, float, num)
    assert ln.Record.filter(feature_int__lt=21).one_or_none() is None
    assert len(ln.Record.filter(feature_int__gt=21)) >= 1
    # int __lt/__gt that would fail with string comparison (42 vs 5, 42 vs 100)
    assert ln.Record.filter(feature_int__lt=5).one_or_none() is None
    assert ln.Record.filter(feature_int__gt=100).one_or_none() is None
    # float/num __lt/__gt (numeric comparison on SQLite via json_extract + CAST)
    assert ln.Record.filter(feature_float__lt=5.0).one() == test_record
    assert ln.Record.filter(feature_float__gt=1.0).one() == test_record
    assert ln.Record.filter(feature_float__gt=10.0).one_or_none() is None
    assert ln.Record.filter(feature_num__lt=5.0).one() == test_record
    assert ln.Record.filter(feature_num__gt=1.0).one() == test_record
    assert ln.Record.filter(feature_num__gt=10.0).one_or_none() is None
    # Date and datetime comparators (ISO strings)
    assert ln.Record.filter(feature_date__lt="2024-01-02").one() == test_record
    assert ln.Record.filter(feature_date__gt="2023-12-31").one() == test_record
    assert ln.Record.filter(feature_date__gt="2024-01-02").one_or_none() is None
    assert (
        ln.Record.filter(feature_datetime__lt="2024-01-01T13:00:00").one()
        == test_record
    )
    assert (
        ln.Record.filter(feature_datetime__gt="2024-01-01T11:00:00").one()
        == test_record
    )
    assert (
        ln.Record.filter(feature_datetime__lt="2024-01-01T11:00:00").one_or_none()
        is None
    )

    # ManyToMany accesors

    assert set(test_record.linked_records.to_list()) == {record_entity1, record_entity2}
    assert test_record.linked_in_records.count() == 0
    assert set(record_entity1.linked_in_records.to_list()) == {test_record}
    assert set(record_entity2.linked_in_records.to_list()) == {test_record}
    assert record_entity1.linked_records.count() == 0
    assert record_entity2.linked_records.count() == 0

    # all empty sheet

    schema = ln.Schema(
        [
            feature_bool,
            feature_str,
            feature_int,
            feature_list_str,
            feature_list_int,
            feature_num,
            feature_float,
            feature_list_float,
            feature_list_num,
            feature_datetime,
            feature_date,
            feature_dict,
            feature_type1,
            feature_type1s,
            feature_ulabel,
            feature_user,
            feature_project,
            feature_cell_line,
            feature_cell_lines,
            feature_cl_ontology_id,
            feature_gene,
            feature_artifact,
            feature_collection,
            feature_run,
        ],
        name="test_schema",
    ).save()
    sheet = ln.Record(name="Sheet", is_type=True, schema=schema).save()
    empty_record = ln.Record(name="empty_record", type=sheet).save()
    df_empty = sheet.to_dataframe()

    assert df_empty["feature_bool"].isnull().all()
    assert df_empty["feature_bool"].dtype.name == "boolean"
    assert df_empty["feature_str"].isnull().all()
    assert df_empty["feature_str"].dtype.name == "string"
    assert df_empty["feature_int"].isnull().all()
    assert df_empty["feature_int"].dtype.name == "Int64"
    assert df_empty["feature_float"].isnull().all()
    assert df_empty["feature_float"].dtype.name == "float64"
    assert df_empty["feature_num"].isnull().all()
    assert df_empty["feature_num"].dtype.name == "float64"
    assert df_empty["feature_list_str"].isnull().all()
    assert df_empty["feature_list_str"].dtype.name == "object"
    assert df_empty["feature_list_int"].isnull().all()
    assert df_empty["feature_list_int"].dtype.name == "object"
    assert df_empty["feature_datetime"].isnull().all()
    assert df_empty["feature_datetime"].dtype.name == "datetime64[ns]"
    assert df_empty["feature_date"].isnull().all()
    assert df_empty["feature_date"].dtype.name == "object"
    assert df_empty["feature_dict"].isnull().all()
    assert df_empty["feature_dict"].dtype.name == "object"
    assert df_empty["feature_type1"].isnull().all()
    assert df_empty["feature_type1"].dtype.name == "category"
    assert df_empty["feature_type1s"].isnull().all()
    assert df_empty["feature_type1s"].dtype.name == "object"
    assert df_empty["feature_ulabel"].isnull().all()
    assert df_empty["feature_ulabel"].dtype.name == "category"
    assert df_empty["feature_user"].isnull().all()
    assert df_empty["feature_user"].dtype.name == "category"
    assert df_empty["feature_project"].isnull().all()
    assert df_empty["feature_project"].dtype.name == "category"
    assert df_empty["feature_cell_line"].isnull().all()
    assert df_empty["feature_cell_line"].dtype.name == "category"
    assert df_empty["feature_cell_lines"].isnull().all()
    assert df_empty["feature_cell_lines"].dtype.name == "object"
    assert df_empty["feature_cl_ontology_id"].isnull().all()
    assert df_empty["feature_cl_ontology_id"].dtype.name == "category"
    assert df_empty["feature_artifact"].isnull().all()
    assert df_empty["feature_artifact"].dtype.name == "category"
    assert df_empty["feature_collection"].isnull().all()
    assert df_empty["feature_collection"].dtype.name == "category"
    assert df_empty["feature_run"].isnull().all()
    assert df_empty["feature_run"].dtype.name == "category"

    # remove empty record from sheet
    empty_record.type = None
    empty_record.save()

    # sheet with values

    test_record.type = sheet
    test_record.save()
    df = sheet.to_dataframe()
    target_result = {
        "feature_bool": True,
        "feature_str": "00810702-0006",  # this string value could be cast to datetime!
        "feature_list_str": ["a", "list", "of", "strings"],
        "feature_int": 42,
        "feature_list_int": [1, 2, 3],
        "feature_float": 3.14,
        "feature_list_float": [2.71, 3.14, 1.61],
        "feature_num": 3.14,
        "feature_list_num": [2.71, 3.14, 1.61],
        "feature_datetime": pd.Timestamp("2024-01-01 12:00:00"),
        "feature_date": date(2024, 1, 1),
        "feature_dict": {"key": "value", "list": [1, 2, 3], "number": 123},
        "feature_type1": "entity1",
        "feature_ulabel": "test-ulabel",
        "feature_user": ln.setup.settings.user.handle,
        "feature_project": "test_project",
        "feature_cell_line": "HEK293",
        "feature_cl_ontology_id": "CVCL_0045",
        "feature_gene": "Tmem276",
        "feature_artifact": "test-artifact",
        "feature_collection": "test-collection",
        "feature_run": run.uid,
        "__lamindb_record_uid__": test_record.uid,
        "__lamindb_record_name__": "test_record",
    }
    result = df.to_dict(orient="records")[0]
    # need to handle categorical lists differently because
    # we don't yet respect ordering
    result_feature_type1s = result.pop("feature_type1s")
    assert set(result_feature_type1s) == {"entity1", "entity2"}
    assert isinstance(result_feature_type1s, list)
    result_feature_cell_lines = result.pop("feature_cell_lines")
    assert set(result_feature_cell_lines) == {"HEK293", "A-549"}
    assert isinstance(result_feature_cell_lines, list)
    assert result == target_result

    # export to artifact to trigger validation -- this will raise many errors if anything is inconsistent

    sheet_as_artifact = sheet.to_artifact()

    # could devise a test for get_values or features.describe()
    # but this is extensively tested elsewhere
    # print(sheet_as_artifact.features.get_values())
    # assert sheet_as_artifact.features.get_values()

    sheet_as_artifact.delete(permanent=True)

    # add the empty record back to the sheet and export again

    empty_record.type = sheet
    empty_record.save()
    df = sheet.to_dataframe()
    sheet_as_artifact = sheet.to_artifact()
    sheet_as_artifact.delete(permanent=True)

    # test passing ISO-format date string for date

    test_record2 = ln.Record(name="test_record").save()
    # we could also test different ways of formatting but don't yet do that
    # in to_dataframe() we enforce ISO format already
    feature_date = ln.Feature.get(name="feature_date")
    feature_date.coerce = True  # have to allow coercion because we're passing a string
    feature_date.save()
    test_values["feature_date"] = "2024-01-02"
    test_record2.features.add_values(test_values)
    test_record2.type = sheet
    test_record2.save()
    test_values["feature_date"] = date(2024, 1, 2)
    assert test_record2.features.get_values() == test_values
    assert test_record.features.get_values() != test_values

    # also test export to artifact again
    sheet_as_artifact = sheet.to_artifact()
    sheet_as_artifact.delete(permanent=True)
    test_record2.delete(permanent=True)
    empty_record.delete(permanent=True)

    # test move a value into the trash

    record_entity1.delete()
    test_values.pop("feature_type1")
    test_values["feature_type1s"] = ["entity2"]
    test_values["feature_date"] = date(2024, 1, 1)
    assert test_record.features.get_values() == test_values

    df = sheet.to_dataframe()
    result = df.to_dict(orient="records")[0]
    result_feature_type1s = result.pop("feature_type1s")
    assert set(result_feature_type1s) == {"entity2"}
    assert isinstance(result_feature_type1s, list)
    result_feature_cell_lines = result.pop("feature_cell_lines")
    assert set(result_feature_cell_lines) == {"HEK293", "A-549"}
    assert isinstance(result_feature_cell_lines, list)
    target_result.pop("feature_type1")
    assert pd.isna(result.pop("feature_type1"))
    assert result == target_result

    record_entity1.restore()
    test_values["feature_type1"] = "entity1"
    test_values["feature_type1s"] = ["entity1", "entity2"]

    # remove values

    test_record.features.remove_values("feature_int")
    test_values.pop("feature_int")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_date")
    test_values.pop("feature_date")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_type1")
    test_values.pop("feature_type1")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_type1s")
    test_values.pop("feature_type1s")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_ulabel")
    test_values.pop("feature_ulabel")
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

    test_record.features.remove_values("feature_collection")
    test_values.pop("feature_collection")
    assert test_record.features.get_values() == test_values

    test_record.features.remove_values("feature_run")
    test_values.pop("feature_run")
    assert test_record.features.get_values() == test_values

    # test passing None has no effect, does not lead to annotation

    sheet.schema = None
    sheet.save()
    schema.delete(permanent=True)

    test_record.features.add_values({"feature_int": None, "feature_type1": None})
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

    # test with list of strings

    schema = ln.Schema([feature_cell_lines], name="test_schema2").save()
    test_form = ln.Record(name="TestForm", is_type=True, schema=schema).save()
    test_record_in_form = ln.Record(name="test_record_in_form", type=test_form).save()
    test_record_in_form.features.add_values({"feature_cell_lines": ["HEK293", "A-549"]})
    test_record_in_form.delete(permanent=True)
    test_form.delete(permanent=True)
    schema.delete(permanent=True)

    # test with list of records (rather than passing strings)

    schema = ln.Schema([feature_cell_lines], name="test_schema2").save()
    test_form = ln.Record(name="TestForm", is_type=True, schema=schema).save()
    test_record_in_form = ln.Record(name="test_record_in_form", type=test_form).save()
    test_record_in_form.features.add_values({"feature_cell_lines": [a549, hek293]})
    test_record_in_form.delete(permanent=True)
    test_form.delete(permanent=True)
    schema.delete(permanent=True)

    # clean up rest
    test_record_id = test_record.id
    assert ln.models.RecordJson.filter(record_id=test_record_id).count() > 0
    test_record.delete(permanent=True)
    # test CASCADE deletion of RecordJson
    assert ln.models.RecordJson.filter(record_id=test_record_id).count() == 0
    sheet.delete(permanent=True)
    feature_str.delete(permanent=True)
    feature_list_str.delete(permanent=True)
    feature_int.delete(permanent=True)
    feature_list_int.delete(permanent=True)
    feature_datetime.delete(permanent=True)
    feature_date.delete(permanent=True)
    feature_type1.delete(permanent=True)
    feature_type1s.delete(permanent=True)
    feature_ulabel.delete(permanent=True)
    feature_user.delete(permanent=True)
    feature_project.delete(permanent=True)
    feature_dict.delete(permanent=True)
    feature_artifact.delete(permanent=True)
    feature_run.delete(permanent=True)
    feature_cell_lines.delete(permanent=True)
    record_entity1.delete(permanent=True)
    record_entity2.delete(permanent=True)
    record_type1.delete(permanent=True)
    test_project.delete(permanent=True)
    feature_cell_line.delete(permanent=True)
    feature_cl_ontology_id.delete(permanent=True)
    feature_collection.delete(permanent=True)
    feature_gene.delete(permanent=True)
    hek293.delete(permanent=True)
    a549.delete(permanent=True)
    tmem276.delete(permanent=True)
    ulabel.delete(permanent=True)
    collection.delete(permanent=True)
    artifact.delete(permanent=True)
    run.delete(permanent=True)
    transform.delete(permanent=True)
    feature_num.delete(permanent=True)


def test_date_and_datetime_corruption():
    feature_datetime = ln.Feature(
        name="feature_datetime", dtype=datetime, coerce=True
    ).save()
    feature_date = ln.Feature(
        name="feature_date", dtype=datetime.date, coerce=True
    ).save()
    schema = ln.Schema(
        [feature_datetime, feature_date], name="test_schema_date_datetime"
    ).save()
    test_sheet = ln.Record(name="TestSheet", is_type=True).save()
    record = ln.Record(name="test_record", type=test_sheet).save()

    # pass values with Z suffix
    test_values = {
        "feature_datetime": "2024-01-01T12:00:00Z",
        "feature_date": "2025-01-17",
    }
    record.features.add_values(test_values)
    date_value = ln.models.RecordJson.get(record=record, feature=feature_date)
    # manually corrupt the value
    date_value.value = "2025-01-17T00:00:00.000Z"
    date_value.save()
    assert record.features.get_values() == {
        "feature_datetime": pd.Timestamp("2024-01-01 12:00:00", tz="UTC"),
        "feature_date": date(2025, 1, 17),
    }
    record.schema = schema
    record.save()

    df = test_sheet.to_dataframe()
    result = df.to_dict(orient="records")[0]
    # because in a dataframe we'll hit pandera and pandera expects naive
    # timestamps, to_dataframe() converts to naive by removing timezone info
    assert result["feature_datetime"] == pd.Timestamp("2024-01-01 12:00:00")
    assert result["feature_date"] == date(2025, 1, 17)

    record.delete(permanent=True)
    test_sheet.delete(permanent=True)
    schema.delete(permanent=True)
    feature_datetime.delete(permanent=True)
    feature_date.delete(permanent=True)


def test_only_list_type_features_and_field_qualifiers():
    # this test is necessary because the logic for adding link tables
    # to the query previously only fired when a non-list cat feature of the same type was present
    feature_cell_lines = ln.Feature(
        name="feature_cell_lines", dtype=list[bt.CellLine]
    ).save()
    feature_list_ontology_id = ln.Feature(
        name="feature_list_ontology_id", dtype=list[bt.Tissue.ontology_id]
    ).save()
    schema = ln.Schema(
        [feature_cell_lines, feature_list_ontology_id], name="test_schema2"
    ).save()
    # create a feature with the same name to test robustness w.r.t. to this
    feature_type = ln.Feature(name="FeatureTypeX", is_type=True).save()
    feature_cell_lines_duplicate = ln.Feature(
        name="feature_cell_lines", dtype=bt.CellLine, type=feature_type
    ).save()

    test_sheet = ln.Record(name="TestSheet", is_type=True, schema=schema).save()
    record = ln.Record(name="test_record", type=test_sheet).save()
    hek293 = bt.CellLine.from_source(name="HEK293").save()
    a549 = bt.CellLine.from_source(name="A-549").save()
    uberon2369 = bt.Tissue.from_source(ontology_id="UBERON:0002369").save()
    uberon5172 = bt.Tissue.from_source(ontology_id="UBERON:0005172").save()

    test_values = {
        "feature_cell_lines": ["HEK293", "A-549"],
        "feature_list_ontology_id": ["UBERON:0002369", "UBERON:0005172"],
    }

    record.features.add_values(test_values)
    assert record.features.get_values() == test_values

    df = test_sheet.to_dataframe()
    result = df.to_dict(orient="records")[0]
    assert isinstance(result["feature_cell_lines"], list)
    assert isinstance(result["feature_list_ontology_id"], list)
    assert set(result["feature_cell_lines"]) == {"HEK293", "A-549"}
    assert set(result["feature_list_ontology_id"]) == {
        "UBERON:0002369",
        "UBERON:0005172",
    }

    # add another record
    record2 = ln.Record(name="test_record2", type=test_sheet).save()
    test_values2 = {
        "feature_cell_lines": ["HEK293"],
        "feature_list_ontology_id": ["UBERON:0005172"],
    }
    record2.features.add_values(test_values2)

    # trigger validation of the case that has two and a single record
    # this tests type casting in list-like values
    artifact = test_sheet.to_artifact()
    assert (
        len(artifact.schemas.first().members) == 2
    )  # this requires top most match filtering during validation

    record.delete(permanent=True)
    record2.delete(permanent=True)
    test_sheet.delete(permanent=True)
    inferred_schema = artifact.schemas.first()
    artifact.delete(permanent=True)
    inferred_schema.delete(permanent=True)
    schema.delete(permanent=True)
    feature_cell_lines.delete(permanent=True)
    feature_cell_lines_duplicate.delete(permanent=True)
    feature_type.delete(permanent=True)
    hek293.delete(permanent=True)
    a549.delete(permanent=True)
    uberon2369.delete(permanent=True)
    uberon5172.delete(permanent=True)
