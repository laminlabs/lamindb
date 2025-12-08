import lamindb as ln
import pandas as pd
import pytest


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


def test_curators_df_nested_cat(nested_cat_df, nested_cat_schema):
    # note that there are two type records both called "Sample" but having different parent types
    # first we create lab_a_type->Sample
    lab_a_type = ln.Record(name="LabA", is_type=True).save()
    sample_a_type = ln.Record(name="Sample", is_type=True, type=lab_a_type).save()
    for name in ["sample1", "sample2", "sample3", "sample4"]:
        ln.Record(name=name, type=sample_a_type).save()

    # then we create lab_b_type->Sample
    lab_b_type = ln.Record(name="LabB", is_type=True).save()
    sample_b_type = ln.Record(name="Sample", is_type=True, type=lab_b_type).save()
    for name in ["sample5", "sample6"]:
        ln.Record(name=name, type=sample_b_type).save()

    # "sample5" is not part of lab_a_type, so it should not be validated
    with pytest.raises(ln.errors.ValidationError):
        curator = ln.curators.DataFrameCurator(nested_cat_df, nested_cat_schema)
        curator.validate()

    assert len(curator.cat._cat_vectors["biosample_name"]._validated) == 4
    assert len(curator.cat._cat_vectors["biosample_name"]._non_validated) == 2

    sample_a_type.records.all().delete(permanent=True)
    sample_b_type.records.all().delete(permanent=True)
    lab_b_type.records.all().delete(permanent=True)
    lab_a_type.records.all().delete(permanent=True)
    lab_a_type.delete(permanet=True)
    lab_b_type.delete(permanent=True)
