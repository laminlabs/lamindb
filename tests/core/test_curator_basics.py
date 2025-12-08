import lamindb as ln
import pandas as pd
import pytest


def test_df_curator_typed_categorical():
    # root level
    sample_root_type = ln.Record(name="Sample", is_type=True).save()
    for name in ["sample1", "sample2"]:
        ln.Record(name=name, type=sample_root_type).save()

    # lab A level
    lab_a_type = ln.Record(name="LabA", is_type=True).save()
    sample_a_type = ln.Record(name="Sample", is_type=True, type=lab_a_type).save()
    for name in ["sample3", "sample4"]:
        ln.Record(name=name, type=sample_a_type).save()

    # lab B level
    lab_b_type = ln.Record(name="LabB", is_type=True).save()
    sample_b_type = ln.Record(name="Sample", is_type=True, type=lab_b_type).save()
    for name in ["sample5", "sample6"]:
        ln.Record(name=name, type=sample_b_type).save()

    feature = ln.Feature(name="biosample_name", dtype=sample_a_type).save()

    df = pd.DataFrame(
        {
            "biosample_name": pd.Categorical(
                ["sample1", "sample2", "sample3", "sample4", "sample5", "sample6"]
            ),
        }
    )

    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert (
        "4 terms not validated in feature 'biosample_name': 'sample1', 'sample2', 'sample5', 'sample6'"
        in error.exconly()
    )
    assert set(curator.cat._cat_vectors["biosample_name"]._validated) == {
        "sample3",
        "sample4",
    }
    assert set(curator.cat._cat_vectors["biosample_name"]._non_validated) == {
        "sample1",
        "sample2",
        "sample5",
        "sample6",
    }

    sample_a_type.records.all().delete(permanent=True)
    sample_b_type.records.all().delete(permanent=True)
    lab_b_type.records.all().delete(permanent=True)
    lab_a_type.records.all().delete(permanent=True)
    lab_a_type.delete(permanet=True)
    lab_b_type.delete(permanent=True)
    feature.delete(permanent=True)
