import bionty as bt
import pandas as pd
import pytest
from lamindb import ULabel
from lamindb.errors import ValidationError
from lamindb.models.feature import (
    parse_dtype,
    parse_filter_expressions,
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
    source.uid = "sourcetestuid"
    source.save()
    return source


@pytest.fixture
def organism():
    organism = bt.Organism(name="test_organism")
    organism.uid = "organismtestuid"
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


# -----------------------------------------------------------------------------
# parsing django filter expressions
# -----------------------------------------------------------------------------


def test_no_registry():
    result = parse_filter_expressions("parent__id=123, category__name=electronics")
    assert result == {"parent__id": "123", "category__name": "electronics"}


def test_relation_filter_with_uid(source):
    result = parse_filter_expressions("source__uid=testuid", bt.Gene)
    assert result == {"source": source}

    source.delete()


def test_relation_filter_with_name(organism):
    result = parse_filter_expressions("organism__name=test_organism", bt.Gene)
    assert result == {"organism": organism}

    organism.delete()


def test_multiple_relation_filters(organism, source):
    result = parse_filter_expressions(
        "organism__name=test_organism, source__uid=sourcetestuid", bt.Gene
    )
    assert result == {"organism": organism, "source": source}

    source.delete()
    organism.delete()


def test_nested_filter(organism):
    result = parse_filter_expressions("organism__name__contains=test_orga", bt.Gene)
    assert result == {"organism": organism}
    organism.delete()


def test_relation_filter_failed_resolution():
    with pytest.raises(bt.Organism.DoesNotExist):
        parse_filter_expressions("organism__name=nonexistent", bt.Gene)


def test_empty_filter():
    result = parse_filter_expressions("", bt.Gene)
    assert result == {}


def test_malformed_filter_no_equals():
    result = parse_filter_expressions("malformed_filter", bt.Gene)
    assert result == {}
