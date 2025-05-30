import bionty
import pandas as pd
import pytest
from lamindb import ULabel
from lamindb.errors import ValidationError
from lamindb.models.feature import parse_dtype, serialize_dtype

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
        "registry": bionty.CellType,
        "field": bionty.CellType.ontology_id,
    }


def test_bionty_perturbations_with_field():
    dtype_str = "cat[bionty.CellType.uid|bionty.CellLine.uid]"
    result = parse_dtype(dtype_str)
    assert len(result) == 2
    assert result[0] == {
        "registry_str": "bionty.CellType",
        "subtype_str": "",
        "field_str": "uid",
        "registry": bionty.CellType,
        "field": bionty.CellType.uid,
    }
    assert result[1] == {
        "registry_str": "bionty.CellLine",
        "subtype_str": "",
        "field_str": "uid",
        "registry": bionty.CellLine,
        "field": bionty.CellLine.uid,
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
    assert serialize_dtype(list[bionty.CellLine]) == "list[cat[bionty.CellLine]]"
