import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb import ULabel
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


def test_seralize_dtypes():
    df = pd.DataFrame(
        {
            "column1": pd.Series([1, 4, 0, 10, 9], dtype="uint"),
        }
    )
    assert df.column1.dtype.name == "uint64"
    assert serialize_dtype(df.column1.dtype) == "int"


# -----------------------------------------------------------------------------
# parsing serialized dtypes
# -----------------------------------------------------------------------------


def test_simple_ulabel_with_subtype_and_field():
    dtype_str = "cat[ULabel[Customer].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "ULabel",
        "subtype_str": "Customer",
        "field_str": "name",
        "registry": ULabel,
        "field": ULabel.name,
    }


def test_multiple_ulabels_with_subtypes_and_fields():
    dtype_str = "cat[ULabel[Customer].name|ULabel[Supplier].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 2
    assert result[0] == {
        "registry_str": "ULabel",
        "subtype_str": "Customer",
        "field_str": "name",
        "registry": ULabel,
        "field": ULabel.name,
    }
    assert result[1] == {
        "registry_str": "ULabel",
        "subtype_str": "Supplier",
        "field_str": "name",
        "registry": ULabel,
        "field": ULabel.name,
    }


def test_bionty_celltype_with_field():
    dtype_str = "cat[bionty.CellType.ontology_id]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "bionty.CellType",
        "subtype_str": "",
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
        "subtype_str": "",
        "field_str": "uid",
        "registry": bt.CellType,
        "field": bt.CellType.uid,
    }
    assert result[1] == {
        "registry_str": "bionty.CellLine",
        "subtype_str": "",
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
    dtype_str = "cat[ULabel[Customer.name"
    with pytest.raises(ValueError) as err:
        parse_dtype(dtype_str)
    assert err.exconly().startswith(
        f"ValueError: dtype is '{dtype_str}' but has to be one of"
    )


def test_simple_registry_without_field():
    dtype_str = "cat[ULabel]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "ULabel",
        "subtype_str": "",
        "field_str": "name",
        "registry": ULabel,
        "field": ULabel.name,
    }


def test_registry_with_subtype_no_field():
    dtype_str = "cat[ULabel[Customer]]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "ULabel",
        "subtype_str": "Customer",
        "field_str": "name",
        "registry": ULabel,
        "field": ULabel.name,
    }


def test_list_of_dtypes():
    dtype_str = "list[cat[ULabel[Customer]]]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "ULabel",
        "subtype_str": "Customer",
        "field_str": "name",
        "registry": ULabel,
        "field": ULabel.name,
        "list": True,
    }
    assert serialize_dtype(list[bt.CellLine]) == "list[cat[bionty.CellLine]]"


def test_nested_cat_dtypes():
    dtype_str = "cat[ULabel[Customer[UScustomer]].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "ULabel",
        "subtype_str": "Customer[UScustomer]",
        "field_str": "name",
        "registry": ULabel,
        "field": ULabel.name,
        "nested_subtypes": ["Customer", "UScustomer"],
    }


def test_registry_with_filter():
    dtype_str = "cat[bionty.Gene.ensembl_gene_id[source__id='abcd']]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {
        "registry_str": "bionty.Gene",
        "subtype_str": "source__id='abcd'",
        "field_str": "ensembl_gene_id",
        "registry": bt.Gene,
        "field": bt.Gene.ensembl_gene_id,
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
        "subtype_str": "source__uid='4a3ejKuf'",
        "field_str": "name",
        "registry": bt.Disease,
        "field": bt.Disease.name,
    }

    feature.delete()


def test_cat_filters_incompatible_with_union_dtypes():
    with pytest.raises(ValidationError) as exc_info:
        ln.Feature(
            name="test_feature",
            dtype="cat[ULabel|bionty.CellType]",
            cat_filters={"source": "test"},
        )
    assert (
        "cat_filters are incompatible with union dtypes: 'cat[ULabel|bionty.CellType]'"
        in str(exc_info.value)
    )


def test_cat_filters_incompatible_with_nested_dtypes():
    with pytest.raises(ValidationError) as exc_info:
        ln.Feature(
            name="test_feature",
            dtype="cat[ULabel[Customer[SubCustomer]]]",
            cat_filters={"source": "test"},
        )
    assert (
        "cat_filters are incompatible with nested dtypes: 'cat[ULabel[Customer[SubCustomer]]]'"
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
    source.delete()


def test_resolve_relation_filter_with_name(organism):
    parsed = {"organism__name": ("organism", "name", "test_organism")}
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"organism": organism}
    organism.delete()


def test_resolve_multiple_relation_filters(organism, source):
    parsed = {
        "organism__name": ("organism", "name", "test_organism"),
        "source__uid": ("source", "uid", "testuid1"),
    }
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"organism": organism, "source": source}
    source.delete()
    organism.delete()


def test_resolve_nested_filter(organism):
    parsed = {"organism__name__contains": ("organism", "name__contains", "test_orga")}
    result = resolve_relation_filters(parsed, bt.Gene)
    assert result == {"organism": organism}
    organism.delete()


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
