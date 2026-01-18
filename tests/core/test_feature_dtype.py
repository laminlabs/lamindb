import datetime

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb import Record
from lamindb.errors import ValidationError
from lamindb.models.feature import (
    dtype_as_object,
    parse_dtype,
    parse_filter_string,
    resolve_relation_filters,
    serialize_dtype,
)


@pytest.fixture
def organism():
    organism = bt.Organism(name="test_organism")
    organism.uid = "testuid2"
    organism.save()
    return organism


# -----------------------------------------------------------------------------
# serializing dtypes
# -----------------------------------------------------------------------------


def test_serialize_basic_dtypes():
    assert serialize_dtype(int) == "int"
    assert serialize_dtype(float) == "float"
    assert serialize_dtype(str) == "str"
    assert serialize_dtype(bool) == "bool"
    assert serialize_dtype(dict) == "dict"
    # assert serialize_dtype(bytes) == "bytes"  # not yet supported
    assert serialize_dtype(datetime.datetime) == "datetime"
    assert serialize_dtype(datetime.date) == "date"


def test_serialize_basic_list_dtypes():
    assert serialize_dtype(list[int]) == "list[int]"
    assert serialize_dtype(list[float]) == "list[float]"
    assert serialize_dtype(list[str]) == "list[str]"
    assert serialize_dtype(list[bool]) == "list[bool]"
    assert serialize_dtype(list[dict]) == "list[dict]"
    assert serialize_dtype(list[datetime.datetime]) == "list[datetime]"
    assert serialize_dtype(list[datetime.date]) == "list[date]"


def test_seralize_pandas_numpy_dtypes():
    series = pd.Series([1, 4, 0, 10, 9], dtype="uint")
    assert series.dtype.name == "uint64"
    assert serialize_dtype(series.dtype) == "int"


def test_serialize_user(ccaplog):
    # correct way through Python object and serialize_dtype()
    feature = ln.Feature(name="user_feat", dtype=ln.User)
    assert feature._dtype_str == "cat[User]"
    # legacy way through parse_dtype()
    feature = ln.Feature(name="user_feat", dtype="cat[User]")
    assert (
        "rather than passing a string 'cat[User]' to dtype, consider passing a Python object"
        in ccaplog.text
    )
    assert feature._dtype_str == "cat[User]"


def test_serialize_record_objects():
    insitute_type = ln.Record(name="InstituteA", is_type=True)
    with pytest.raises(ln.errors.InvalidArgument) as error:
        serialize_dtype(insitute_type)
    assert (
        f"Cannot serialize unsaved objects. Save {insitute_type} via `.save()`."
        in error.exconly()
    )
    insitute_type.save()
    lab_type = ln.Record(name="LabB", type=insitute_type, is_type=True).save()
    sample_type = ln.Record(name="Sample", type=lab_type, is_type=True).save()
    # New UID-based format: cat[Record[uid]] instead of cat[Record[Parent[Child]]]
    serialized_str = f"cat[Record[{sample_type.uid}]]"
    feature = ln.Feature(name="sample_feature", dtype=sample_type).save()
    assert feature._dtype_str == serialized_str
    assert feature.dtype == "cat[Record[InstituteA[LabB[Sample]]]]"
    feature.delete(permanent=True)
    assert serialize_dtype(sample_type) == serialized_str
    with pytest.raises(ln.errors.IntegrityError) as error:
        parse_dtype("cat[Record[Sample]]", check_exists=True, old_format=True)
    assert (
        "No Record type found matching subtypes ['Sample'] for field `.name`"
        in error.exconly()
    )
    sample = ln.Record(name="sample").save()
    with pytest.raises(ln.errors.InvalidArgument) as error:
        parse_dtype(f"cat[Record[{sample.uid}]]", check_exists=True)
    assert (
        f"The resolved Record 'sample' (uid='{sample.uid}') is not a type: is_type is False."
        in error.exconly()
    )
    with pytest.raises(ln.errors.InvalidArgument) as error:
        serialize_dtype(sample)
    assert (
        "Cannot serialize non-type Record 'sample'. Only types (is_type=True) are allowed in dtypes."
        in error.exconly()
    )
    sample_type.delete(permanent=True)
    lab_type.delete(permanent=True)
    insitute_type.delete(permanent=True)
    sample.delete(permanent=True)


def test_serialize_union_of_registries():
    serialized_str = "cat[Record|bionty.Gene]"
    assert serialize_dtype([ln.Record, bt.Gene]) == serialized_str
    serialized_str = "cat[bionty.CellType|bionty.CellLine]"
    assert serialize_dtype([bt.CellType, bt.CellLine]) == serialized_str


def test_serialize_with_field_information():
    serialized_str = "cat[bionty.Gene.ensembl_gene_id]"
    assert serialize_dtype(bt.Gene.ensembl_gene_id) == serialized_str
    serialized_str = "cat[bionty.CellType.uid|bionty.CellLine.uid]"
    assert serialize_dtype([bt.CellType.uid, bt.CellLine.uid]) == serialized_str


# -----------------------------------------------------------------------------
# parsing serialized dtypes
# -----------------------------------------------------------------------------


def test_simple_record_with_subtype_and_field():
    # Create a Record type to get its UID
    customer_type = ln.Record(name="Customer", is_type=True).save()
    dtype_str = f"cat[Record[{customer_type.uid}].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "record_uid": customer_type.uid,
    }
    customer_type.delete(permanent=True)


def test_multiple_records_with_subtypes_and_fields():
    # Create Record types to get their UIDs
    customer_type = ln.Record(name="Customer", is_type=True).save()
    supplier_type = ln.Record(name="Supplier", is_type=True).save()
    dtype_str = (
        f"cat[Record[{customer_type.uid}].name|Record[{supplier_type.uid}].name]"
    )
    result = parse_dtype(dtype_str)
    assert len(result) == 2
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "record_uid": customer_type.uid,
    }
    assert result[1] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "record_uid": supplier_type.uid,
    }
    customer_type.delete(permanent=True)
    supplier_type.delete(permanent=True)


def test_bionty_celltype_with_field():
    dtype_str = "cat[bionty.CellType.ontology_id]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "bionty.CellType",
        "filter_str": "",
        "field_str": "ontology_id",
        "registry": bt.CellType,
        "field": bt.CellType.ontology_id,
    }


def test_bionty_perturbations_with_field():
    dtype_str = "cat[bionty.CellType.uid|bionty.CellLine.uid]"
    result = parse_dtype(dtype_str)
    assert len(result) == 2
    assert result[0] == {
        "registry_str": "bionty.CellType",
        "filter_str": "",
        "field_str": "uid",
        "registry": bt.CellType,
        "field": bt.CellType.uid,
    }
    assert result[1] == {
        "registry_str": "bionty.CellLine",
        "filter_str": "",
        "field_str": "uid",
        "registry": bt.CellLine,
        "field": bt.CellLine.uid,
    }


def test_invalid_registry():
    dtype_str = "cat[InvalidRegistry.field]"
    with pytest.raises(ValidationError) as exc_info:
        parse_dtype(dtype_str)
    assert "invalid dtype" in str(exc_info.value)


def test_empty_category():
    dtype_str = "cat[]"
    result = parse_dtype(dtype_str)
    assert result == []


def test_malformed_categorical():
    dtype_str = "cat ? str"
    with pytest.raises(ValueError) as err:
        parse_dtype(dtype_str)
    assert err.exconly().startswith(
        f"ValueError: dtype is '{dtype_str}' but has to be one of"
    )
    dtype_str = "cat[Record[Customer.name"
    with pytest.raises(ValueError) as err:
        parse_dtype(dtype_str)
    assert err.exconly().startswith(
        f"ValueError: dtype is '{dtype_str}' but has to be one of"
    )


def test_simple_registry_without_field():
    dtype_str = "cat[Record]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
    }


def test_registry_with_subtype_no_field():
    # Create a Record type to get its UID
    customer_type = ln.Record(name="Customer", is_type=True).save()
    dtype_str = f"cat[Record[{customer_type.uid}]]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "record_uid": customer_type.uid,
    }
    customer_type.delete(permanent=True)


def test_list_of_dtypes():
    # Create a Record type to get its UID
    customer_type = ln.Record(name="Customer", is_type=True).save()
    dtype_str = f"list[cat[Record[{customer_type.uid}]]]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "record_uid": customer_type.uid,
        "list": True,
    }
    assert serialize_dtype(list[bt.CellLine]) == "list[cat[bionty.CellLine]]"
    customer_type.delete(permanent=True)


def test_registry_with_filter():
    dtype_str = "cat[bionty.Gene.ensembl_gene_id[source__id='abcd']]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "bionty.Gene",
        "filter_str": "source__id='abcd'",
        "field_str": "ensembl_gene_id",
        "registry": bt.Gene,
        "field": bt.Gene.ensembl_gene_id,
    }


def test_nested_cat_dtypes():
    # Create Record types - the deepest type is UScustomer
    customer_type = ln.Record(name="Customer", is_type=True).save()
    uscustomer_type = ln.Record(
        name="UScustomer", type=customer_type, is_type=True
    ).save()
    dtype_str = f"cat[Record[{uscustomer_type.uid}].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "record_uid": uscustomer_type.uid,
    }
    uscustomer_type.delete(permanent=True)
    customer_type.delete(permanent=True)


def test_nested_cat_with_filter():
    # Create Record types - the deepest type is UScustomer
    # Note: filters in bracket content are not currently supported in UID format
    # This test may need adjustment based on how filters are handled
    customer_type = ln.Record(name="Customer", is_type=True).save()
    uscustomer_type = ln.Record(
        name="UScustomer", type=customer_type, is_type=True
    ).save()
    dtype_str = f"cat[Record[{uscustomer_type.uid}].description]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "description",
        "registry": Record,
        "field": Record.description,
        "record_uid": uscustomer_type.uid,
    }
    uscustomer_type.delete(permanent=True)
    customer_type.delete(permanent=True)


# -----------------------------------------------------------------------------
# parsing django filter expressions
# -----------------------------------------------------------------------------


def test_feature_dtype():
    feature = ln.Feature(
        name="disease",
        dtype=bt.Disease,
        cat_filters={
            "source__uid": "4a3ejKuf"
        },  # uid corresponds to disease_ontology_old.uid
    ).save()

    result = parse_dtype(feature._dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "bionty.Disease",
        "filter_str": "source__uid='4a3ejKuf'",
        "field_str": "name",
        "registry": bt.Disease,
        "field": bt.Disease.name,
    }

    feature.delete(permanent=True)


def test_cat_filters_incompatible_with_union_dtypes():
    with pytest.raises(ValidationError) as exc_info:
        ln.Feature(
            name="test_feature",
            dtype="cat[Record|bionty.CellType]",
            cat_filters={"source": "test"},
        )
    assert (
        "cat_filters are incompatible with union dtypes: 'cat[Record|bionty.CellType]'"
        in str(exc_info.value)
    )


def test_cat_filters_incompatible_with_nested_dtypes():
    record = ln.Record(name="Customer", is_type=True).save()
    with pytest.raises(ValidationError) as exc_info:
        ln.Feature(
            name="test_feature",
            dtype=record,
            cat_filters={"source": "test"},
        )
    assert (
        f"cat_filters are incompatible with nested dtypes: 'cat[Record[{record.uid}]]'"
        in str(exc_info.value)
    )
    record.delete(permanent=True)


def test_parse_filter_string_basic():
    result = parse_filter_string("parent__id=123, category__name=electronics")
    expected = {
        "parent__id": ("parent", "id", "123"),
        "category__name": ("category", "name", "electronics"),
    }
    assert result == expected


def test_parse_filter_string_direct_fields():
    result = parse_filter_string("name=test, status=active")
    expected = {"name": ("name", None, "test"), "status": ("status", None, "active")}
    assert result == expected


def test_parse_filter_string_empty():
    with pytest.raises(ValueError) as e:
        parse_filter_string("")
        assert "missing '=' sign" in str(e)


def test_parse_filter_string_malformed():
    with pytest.raises(ValueError) as e:
        parse_filter_string("malformed_filter")
        assert "missing '=' sign" in str(e)


def test_parse_filter_string_missing_key():
    with pytest.raises(ValueError) as e:
        parse_filter_string("=someval")
        assert "empty key" in str(e)


def test_parse_filter_string_missing_value():
    with pytest.raises(ValueError) as e:
        parse_filter_string("somekey=")
        assert "empty val" in str(e)


def test_resolve_direct_fields():
    parsed = {"name": ("name", None, "test"), "status": ("status", None, "active")}
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"name": "test", "status": "active"}


def test_resolve_relation_filter_with_uid():
    source = bt.Source(
        name="test_name",
        description="test_description",
        organism="human",
        entity="bionty.Gene",
        version="2026-01-01",
    )
    source.uid = "testuid1"
    source.save()
    parsed = {"source__uid": ("source", "uid", "testuid1")}
    result = resolve_relation_filters(parsed, bt.Gene)
    print(result)
    assert result == {"source": source}
    source.delete(permanent=True)


def test_resolve_relation_filter_with_name(organism):
    parsed = {"organism__name": ("organism", "name", "test_organism")}
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"organism": organism}
    organism.delete(permanent=True)


def test_resolve_multiple_relation_filters(organism):
    source = bt.Source(
        name="test_name",
        description="test_description",
        organism="human",
        entity="bionty.Gene",
        version="2026-01-01",
    )
    source.uid = "testuid1"
    source.save()
    parsed = {
        "organism__name": ("organism", "name", "test_organism"),
        "source__uid": ("source", "uid", "testuid1"),
    }
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"organism": organism, "source": source}
    source.delete(permanent=True)
    organism.delete(permanent=True)


def test_resolve_nested_filter(organism):
    parsed = {"organism__name__contains": ("organism", "name__contains", "test_orga")}
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"organism": organism}
    organism.delete(permanent=True)


def test_resolve_relation_filter_failed_resolution():
    parsed = {"organism__name": ("organism", "name", "nonexistent")}
    with pytest.raises(bt.Organism.DoesNotExist):
        resolve_relation_filters(parsed, bt.Gene)


def test_resolve_relation_filter_duplicate():
    parsed = {
        "source__uid": ("source", "uid", "testuid1"),
        "source__name": ("source", "name", "test_name"),
    }
    with pytest.raises(bt.Source.DoesNotExist):
        resolve_relation_filters(parsed, bt.Gene)


# -----------------------------------------------------------------------------
# backward compatibility for old format strings
# -----------------------------------------------------------------------------


def test_convert_old_format_ulabel_string():
    """Test converting old format ULabel string to object."""
    # Create a ULabel type
    perturbation = ln.ULabel(name="Perturbation", is_type=True).save()

    # Convert old format string
    dtype_obj = dtype_as_object("cat[ULabel[Perturbation]]", old_format=True)

    # Should return the ULabel object
    assert dtype_obj == perturbation
    assert hasattr(dtype_obj, "uid")

    # Clean up
    perturbation.delete(permanent=True)


def test_convert_old_format_record_string():
    """Test converting old format Record string to object."""
    # Create a Record type
    sample_type = ln.Record(name="Sample", is_type=True).save()

    # Convert old format string
    dtype_obj = dtype_as_object("cat[Record[Sample]]", old_format=True)

    # Should return the Record object
    assert dtype_obj == sample_type
    assert hasattr(dtype_obj, "uid")

    # Clean up
    sample_type.delete(permanent=True)


def test_convert_old_format_nested_record_string():
    """Test converting old format nested Record string to object."""
    # Create nested Record types
    lab_type = ln.Record(name="LabA", is_type=True).save()
    experiment_type = ln.Record(name="Experiment", type=lab_type, is_type=True).save()

    # Convert old format string
    dtype_obj = dtype_as_object("cat[Record[LabA[Experiment]]]", old_format=True)

    # Should return the nested Record object
    assert dtype_obj == experiment_type
    assert hasattr(dtype_obj, "uid")

    # Clean up
    experiment_type.delete(permanent=True)
    lab_type.delete(permanent=True)


def test_convert_old_format_list_string():
    """Test converting old format list string to object."""
    # Create a ULabel type
    perturbation = ln.ULabel(name="Perturbation", is_type=True).save()

    # Convert old format string with list wrapper
    dtype_obj = dtype_as_object("list[cat[ULabel[Perturbation]]]", old_format=True)

    # Should return list[ULabel] type
    assert hasattr(dtype_obj, "__origin__")
    assert dtype_obj.__origin__ is list
    # Get the inner type
    from typing import get_args

    inner_type = get_args(dtype_obj)[0]
    assert inner_type == perturbation

    # Clean up
    perturbation.delete(permanent=True)


def test_feature_constructor_with_old_format_string(ccaplog):
    """Test Feature constructor with old format string raises deprecation warning."""
    # Create a ULabel type
    perturbation = ln.ULabel(name="Perturbation", is_type=True).save()

    # Create feature with old format string
    feature = ln.Feature(name="perturbation", dtype="cat[ULabel[Perturbation]]")
    assert (
        "rather than passing a string 'cat[ULabel[Perturbation]]' to dtype, consider passing a Python object"
        in ccaplog.text
    )

    # Should have converted to UID format
    assert feature._dtype_str is not None
    assert "ULabel[" in feature._dtype_str
    # Should contain UID, not name
    assert "Perturbation" not in feature._dtype_str
    assert perturbation.uid in feature._dtype_str

    # Clean up
    perturbation.delete(permanent=True)


def test_feature_constructor_with_old_format_nested_string(ccaplog):
    """Test Feature constructor with old format nested string."""
    # Create nested Record types
    lab_type = ln.Record(name="LabA", is_type=True).save()
    experiment_type = ln.Record(name="Experiment", type=lab_type, is_type=True).save()

    # Create feature with old format nested string
    feature = ln.Feature(name="experiment", dtype="cat[Record[LabA[Experiment]]]")
    assert (
        "rather than passing a string 'cat[Record[LabA[Experiment]]]' to dtype, consider passing a Python object"
        in ccaplog.text
    )

    # Should have converted to UID format
    assert feature._dtype_str is not None
    assert "Record[" in feature._dtype_str
    # Should contain UID, not names
    assert "LabA" not in feature._dtype_str
    assert "Experiment" not in feature._dtype_str
    assert experiment_type.uid in feature._dtype_str

    # Clean up
    experiment_type.delete(permanent=True)
    lab_type.delete(permanent=True)


def test_bare_cat_dtype_backward_compatibility():
    """Test that bare 'cat' dtype is accepted for backward compatibility."""
    # Test parse_dtype accepts "cat" and returns empty list
    result = parse_dtype("cat")
    assert result == []

    # Test Feature constructor with bare "cat" dtype issues deprecation warning
    with pytest.warns(DeprecationWarning, match="dtype `cat` is deprecated"):
        feature = ln.Feature(name="test_bare_cat", dtype="cat")
    assert feature._dtype_str == "cat"


def test_migrate_dtype_to_uid_format():
    """Test migrate_dtype_to_uid_format() function for migration."""
    from django.db import connection
    from lamindb.models.feature import migrate_dtype_to_uid_format

    # Create Record types for testing
    lab_type = ln.Record(name="LabA", is_type=True).save()
    experiment_type = ln.Record(name="Experiment", type=lab_type, is_type=True).save()
    perturbation = ln.ULabel(name="Perturbation", is_type=True).save()

    # Create features with old format strings in _dtype_str
    feature1 = ln.Feature(name="test_record_old_format", dtype="str").save()
    feature2 = ln.Feature(name="test_ulabel_old_format", dtype="str").save()
    feature3 = ln.Feature(name="test_list_record_old_format", dtype="str").save()
    feature4 = ln.Feature(name="test_list_ulabel_old_format", dtype="str").save()

    # Manually set old format strings using raw SQL
    old_format_record = "cat[Record[LabA[Experiment]]]"
    old_format_ulabel = "cat[ULabel[Perturbation]]"
    old_format_list_record = "list[cat[Record[LabA[Experiment]]]]"
    old_format_list_ulabel = "list[cat[ULabel[Perturbation]]]"

    with connection.cursor() as cursor:
        cursor.execute(
            "UPDATE lamindb_feature SET _dtype_str = %s WHERE id = %s",
            [old_format_record, feature1.id],
        )
        cursor.execute(
            "UPDATE lamindb_feature SET _dtype_str = %s WHERE id = %s",
            [old_format_ulabel, feature2.id],
        )
        cursor.execute(
            "UPDATE lamindb_feature SET _dtype_str = %s WHERE id = %s",
            [old_format_list_record, feature3.id],
        )
        cursor.execute(
            "UPDATE lamindb_feature SET _dtype_str = %s WHERE id = %s",
            [old_format_list_ulabel, feature4.id],
        )

    # Refresh features from database
    feature1.refresh_from_db()
    feature2.refresh_from_db()
    feature3.refresh_from_db()
    feature4.refresh_from_db()

    # Verify old format is present
    assert feature1._dtype_str == old_format_record
    assert feature2._dtype_str == old_format_ulabel
    assert feature3._dtype_str == old_format_list_record
    assert feature4._dtype_str == old_format_list_ulabel

    # Run migration function
    migrate_dtype_to_uid_format(connection, input_field="_dtype_str")

    # Refresh features from database
    feature1.refresh_from_db()
    feature2.refresh_from_db()
    feature3.refresh_from_db()
    feature4.refresh_from_db()

    # Verify conversion to UID format
    assert feature1._dtype_str == f"cat[Record[{experiment_type.uid}]]"
    assert feature2._dtype_str == f"cat[ULabel[{perturbation.uid}]]"
    assert feature3._dtype_str == f"list[cat[Record[{experiment_type.uid}]]]"
    assert feature4._dtype_str == f"list[cat[ULabel[{perturbation.uid}]]]"

    # Verify old names are not in the converted strings
    assert "LabA" not in feature1._dtype_str
    assert "Experiment" not in feature1._dtype_str
    assert "Perturbation" not in feature2._dtype_str
    assert "LabA" not in feature3._dtype_str
    assert "Experiment" not in feature3._dtype_str
    assert "Perturbation" not in feature4._dtype_str

    # Verify UIDs are present
    assert experiment_type.uid in feature1._dtype_str
    assert perturbation.uid in feature2._dtype_str
    assert experiment_type.uid in feature3._dtype_str
    assert perturbation.uid in feature4._dtype_str

    # Clean up
    feature1.delete(permanent=True)
    feature2.delete(permanent=True)
    feature3.delete(permanent=True)
    feature4.delete(permanent=True)
    experiment_type.delete(permanent=True)
    lab_type.delete(permanent=True)
    perturbation.delete(permanent=True)
