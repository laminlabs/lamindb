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
def cat_df():
    return pd.DataFrame(
        {
            "sample_id": [["sample1", "sample2"], ["sample2"], ["sample3"]],
            "dose": [[1.2, 2.3], [1.2], [2.3]],
            "cell_type": [["B cell", "T cell"], ["B cell"], ["T cell"]],
            "tissue": ["blood", "blood", "lung"],
        }
    )


@pytest.fixture(scope="module")
def nested_cat_df():
    return pd.DataFrame(
        {
            "biosample_id": ["S1", "S2", "S3", "S4", "S5", "S6"],
            "biosample_name": [
                "sample1",
                "sample2",
                "sample3",
                "sample4",
                "sample5",
                "sample6",
            ],
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

    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    bt.Tissue.filter().delete(permanent=True)


@pytest.fixture(scope="module")
def nested_cat_schema():
    schema = ln.Schema(
        name="schema_with_nested_cat",
        features=[
            ln.Feature(name="biosample_id", dtype=str).save(),
            ln.Feature(name="biosample_name", dtype="cat[Record[LabA[Sample]]]").save(),
        ],
        coerce_dtype=True,
    ).save()

    yield schema

    ln.Schema.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)
    ln.Record.filter().update(type=None)
    ln.Record.filter().delete(permanent=True)


def test_curator_df_multivalue(lists_df, lists_schema, cat_df):
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

    # test with cat_df which has a non-list tissue
    curator = ln.curators.DataFrameCurator(cat_df, lists_schema)
    with pytest.raises(ValidationError):
        curator.validate()


def test_curators_df_nested_cat(nested_cat_df, nested_cat_schema):
    # note that there are two type records both called "Sample" but having different parent types
    # first we create LabA->Sample
    LabA = ln.Record(name="LabA", is_type=True).save()
    Sample = ln.Record(name="Sample", is_type=True, type=LabA).save()
    for name in ["sample1", "sample2", "sample3", "sample4"]:
        ln.Record(name=name, type=Sample).save()

    # then we create LabB->Sample
    LabB = ln.Record(name="LabB", is_type=True).save()
    Sample = ln.Record(name="Sample", is_type=True, type=LabB).save()
    for name in ["sample5", "sample6"]:
        ln.Record(name=name, type=Sample).save()

    # "sample5" is not part of LabA, so it should not be validated
    with pytest.raises(ValidationError):
        curator = ln.curators.DataFrameCurator(nested_cat_df, nested_cat_schema)
        curator.validate()

    assert len(curator.cat._cat_vectors["biosample_name"]._validated) == 4
    assert len(curator.cat._cat_vectors["biosample_name"]._non_validated) == 2


def test_curators_list_feature_nullable_empty_list():
    """Test that a list feature that is nullable can accept empty lists."""
    feature_list = ln.Feature(
        name="list_tissue", dtype=list[bt.Tissue.ontology_id], nullable=True
    ).save()
    feature_int = ln.Feature(name="feature int", dtype=int, nullable=True).save()
    schema = ln.Schema(
        name="test_list_feature_schema",
        features=[feature_list, feature_int],
        coerce_dtype=True,
    ).save()

    df = pd.DataFrame({"list_tissue": [], "feature int": []})
    ln.curators.DataFrameCurator(df, schema).validate()

    # clean up
    schema.delete(permanent=True)
    feature_list.delete(permanent=True)
    feature_int.delete(permanent=True)
