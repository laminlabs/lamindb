import pytest
from lamindb._feature import parse_dtype
from lamindb.core.exceptions import ValidationError


# Test cases
def test_simple_ulabel_with_subtype_and_field():
    dtype_str = "cat[ULabel[Customer].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {"registry": "ULabel", "subtype": "Customer", "field": "name"}


def test_multiple_ulabels_with_subtypes_and_fields():
    dtype_str = "cat[ULabel[Customer].name|ULabel[Supplier].name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 2
    assert result[0] == {"registry": "ULabel", "subtype": "Customer", "field": "name"}
    assert result[1] == {"registry": "ULabel", "subtype": "Supplier", "field": "name"}


def test_bionty_celltype_with_field():
    dtype_str = "cat[bionty.CellType.name]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {"registry": "bionty.CellType", "subtype": "", "field": "name"}


def test_bionty_perturbations_with_field():
    dtype_str = "cat[bionty.CellType.uid|bionty.CellLine.uid]"
    result = parse_dtype(dtype_str)
    assert len(result) == 2
    assert result[0] == {"registry": "bionty.CellType", "subtype": "", "field": "uid"}
    assert result[1] == {"registry": "bionty.CellLine", "subtype": "", "field": "uid"}


def test_invalid_registry():
    dtype_str = "cat[InvalidRegistry.field]"
    with pytest.raises(ValidationError) as exc_info:
        parse_dtype(dtype_str)
    assert "invalid dtype" in str(exc_info.value)


def test_empty_category():
    dtype_str = "cat[]"
    result = parse_dtype(dtype_str)
    assert result == []


def test_malformed_category_missing_bracket():
    dtype_str = "cat[ULabel[Customer.name"
    with pytest.raises(AssertionError):
        parse_dtype(dtype_str)


def test_simple_registry_without_field():
    dtype_str = "cat[ULabel]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {"registry": "ULabel", "subtype": "", "field": ""}


def test_registry_with_subtype_no_field():
    dtype_str = "cat[ULabel[Customer]]"
    result = parse_dtype(dtype_str)
    assert len(result) == 1
    assert result[0] == {"registry": "ULabel", "subtype": "Customer", "field": ""}


def test_not_category_dtype():
    dtype_str = "string"
    result = parse_dtype(dtype_str)
    assert result == []
