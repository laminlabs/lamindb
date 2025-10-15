import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.core.exceptions import ValidationError


@pytest.fixture(scope="module")
def lists_df():
    return pd.DataFrame(
        {
            "sample_id": [["sample1", "sample2"], ["sample2"], ["sample3"]],
            "dose": [[1.2, 2.3], [1.2], [2.3]],
            "cell_type": [["B cell", "T cell"], ["B cell"], ["T cell"]],
            "tissue": [["blood", "pulmo"], ["blood"], ["lung"]],
        }
    )


@pytest.fixture(scope="module")
def nested_cat_df():
    return pd.DataFrame(
        {
            "customer_id": ["C1", "C2", "C3", "C4", "C5", "C6"],
            "us_customer_name": ["Alice", "Bob", "Charlie", "David", "Eve", "Frank"],
        }
    )


@pytest.fixture(scope="module")
def lists_schema():
    schema = ln.Schema(
        name="lists schema cat",
        features=[
            ln.Feature(name="sample_id", dtype=list[str]).save(),
            ln.Feature(name="dose", dtype=list[float]).save(),
            ln.Feature(name="cell_type", dtype=list[str]).save(),
            ln.Feature(name="tissue", dtype=list[bt.Tissue]).save(),
        ],
    ).save()

    yield schema

    schema.delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    bt.Tissue.filter().delete(permanent=True)


@pytest.fixture(scope="module")
def nested_cat_schema():
    schema = ln.Schema(
        name="schema_with_nested_cat",
        features=[
            ln.Feature(name="customer_id", dtype=str).save(),
            ln.Feature(
                name="us_customer_name", dtype="cat[Record[UScustomer[Customer]]]"
            ).save(),
        ],
        coerce_dtype=True,
    ).save()

    yield schema

    schema.delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    ln.Record.filter().update(type=None)
    ln.Record.filter().delete(permanent=True)


def test_curator_df_multivalue(lists_df, lists_schema):
    curator = ln.curators.DataFrameCurator(lists_df, lists_schema)
    with pytest.raises(ValidationError):
        curator.validate()
    assert curator.cat._cat_vectors.keys() == {"columns", "tissue"}
    assert curator.cat._cat_vectors["tissue"]._validated == ["blood", "lung"]
    assert curator.cat._cat_vectors["tissue"]._non_validated == ["pulmo"]
    assert curator.cat._cat_vectors["tissue"]._synonyms == {"pulmo": "lung"}

    curator.cat.standardize("tissue")
    assert curator.cat._cat_vectors["tissue"]._non_validated == []
    assert lists_df["tissue"].tolist() == [["blood", "lung"], ["blood"], ["lung"]]

    assert curator.validate() is None


def test_curators_df_nested_cat(nested_cat_df, nested_cat_schema):
    # Create hierarchical Record types and some Records
    UScustomer = ln.Record(name="UScustomer", is_type=True).save()
    Customer = ln.Record(name="Customer", is_type=True, type=UScustomer).save()
    for name in ["Alice", "Bob", "Charlie", "David"]:
        ln.Record(name=name, type=Customer).save()
    EUcustomer = ln.Record(name="EUcustomer", is_type=True).save()
    Customer = ln.Record(
        name="Customer", is_type=True, type=EUcustomer, _skip_validation=True
    ).save()
    for name in ["Eve", "Frank"]:
        ln.Record(name=name, type=Customer).save()

    # "test_customer" is not a UScustomer, so it should not be validated
    with pytest.raises(ValidationError):
        curator = ln.curators.DataFrameCurator(nested_cat_df, nested_cat_schema)
        curator.validate()

    assert len(curator.cat._cat_vectors["us_customer_name"]._validated) == 4
    assert len(curator.cat._cat_vectors["us_customer_name"]._non_validated) == 2
