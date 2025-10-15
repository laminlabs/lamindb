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
def source():
    source = bt.Source(
        name="test_name",
        description="test_description",
        organism="human",
        entity="bionty.Gene",
        version="2026-01-01",
    )
    source.uid = "testuid1"
    source.save()
    return source


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


def test_serialize_user():
    feature = ln.Feature(
        name="user_feat", dtype="cat[User]"
    )  # calls parse_dtype() to verify
    feature = ln.Feature(name="user_feat", dtype=ln.User)  # calls serialize_dtype()
    assert feature.dtype == "cat[User]"


def test_serialize_record_objects():
    record_type_ist1 = ln.Record(name="Institute1", is_type=True).save()
    record_type_dpt1 = ln.Record(
        name="Department1", type=record_type_ist1, is_type=True
    ).save()
    record_type_lab = ln.Record(
        name="Instrument", type=record_type_dpt1, is_type=True
    ).save()
    serialized_str = "cat[Record[Institute1[Department1[Instrument]]]]"
    assert serialize_dtype(record_type_lab) == serialized_str
    record_type_lab.delete(permanent=True)
    record_type_dpt1.delete(permanent=True)
    record_type_ist1.delete(permanent=True)


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


def test_serialize_with_additional_filters():
    pass
    # see parse_dtype
    # see parse_dtype


# -----------------------------------------------------------------------------
# parsing serialized dtypes
# -----------------------------------------------------------------------------


def test_simple_record_with_subtype_and_field():
    dtype_str = "cat[Record[Customer].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "subtypes_list": ["Customer"],
    }


def test_multiple_records_with_subtypes_and_fields():
    dtype_str = "cat[Record[Customer].name|Record[Supplier].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 2
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "subtypes_list": ["Customer"],
    }
    assert result[1] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "subtypes_list": ["Supplier"],
    }


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
        "subtypes_list": [],
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
        "subtypes_list": [],
    }
    assert result[1] == {
        "registry_str": "bionty.CellLine",
        "filter_str": "",
        "field_str": "uid",
        "registry": bt.CellLine,
        "field": bt.CellLine.uid,
        "subtypes_list": [],
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
        "subtypes_list": [],
    }


def test_registry_with_subtype_no_field():
    dtype_str = "cat[Record[Customer]]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "subtypes_list": ["Customer"],
    }


def test_list_of_dtypes():
    dtype_str = "list[cat[Record[Customer]]]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "subtypes_list": ["Customer"],
        "list": True,
    }
    assert serialize_dtype(list[bt.CellLine]) == "list[cat[bionty.CellLine]]"


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
        "subtypes_list": [],
    }


def test_nested_cat_dtypes():
    dtype_str = "cat[Record[Customer[UScustomer]].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "",
        "field_str": "name",
        "registry": Record,
        "field": Record.name,
        "subtypes_list": ["Customer", "UScustomer"],
    }


def test_nested_cat_with_filter():
    dtype_str = "cat[Record[Customer[UScustomer[region='US']]].description]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "Record",
        "filter_str": "region='US'",
        "field_str": "description",
        "registry": Record,
        "field": Record.description,
        "subtypes_list": ["Customer", "UScustomer"],
    }


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

    result = parse_dtype(feature.dtype)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "bionty.Disease",
        "filter_str": "source__uid='4a3ejKuf'",
        "field_str": "name",
        "registry": bt.Disease,
        "field": bt.Disease.name,
        "subtypes_list": [],
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
    with pytest.raises(ValidationError) as exc_info:
        ln.Feature(
            name="test_feature",
            dtype="cat[Record[Customer[SubCustomer]]]",
            cat_filters={"source": "test"},
        )
    assert (
        "cat_filters are incompatible with nested dtypes: 'cat[Record[Customer[SubCustomer]]]'"
        in str(exc_info.value)
    )


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


def test_resolve_relation_filter_with_uid(source):
    parsed = {"source__uid": ("source", "uid", "testuid1")}
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"source": source}
    source.delete(permanent=True)


def test_resolve_relation_filter_with_name(organism):
    parsed = {"organism__name": ("organism", "name", "test_organism")}
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"organism": organism}
    organism.delete(permanent=True)


def test_resolve_multiple_relation_filters(organism, source):
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
