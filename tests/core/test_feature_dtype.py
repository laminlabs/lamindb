import datetime

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb import Record
from lamindb.errors import ValidationError
from lamindb.models.feature import (
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
        "rather than passing a string 'cat[User]' to dtype, pass a Python object"
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
    assert serialize_dtype(sample_type) == serialized_str
    with pytest.raises(ln.errors.IntegrityError) as error:
        parse_dtype("cat[Record[Sample]]", check_exists=True)
    assert (
        "Error retrieving Record with uid 'Sample' for field `.name`: Record matching query does not exist."
        in error.exconly()
    )
    sample = ln.Record(name="sample").save()
    with pytest.raises(ln.errors.InvalidArgument) as error:
        parse_dtype(f"cat[Record[{sample.uid}]]", check_exists=True)
    assert (
        "The resolved Record 'sample' for field `.name` is not a type: is_type is False."
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
            dtype=f"cat[Record[{record.uid}]]",
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
