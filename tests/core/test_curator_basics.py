import lamindb as ln
import pandas as pd
import pytest


def test_df_curator_typed_categorical():
    # root level
    sample_root_type = ln.Record(name="Sample", is_type=True).save()
    for name in ["s1", "s2"]:
        ln.Record(name=name, type=sample_root_type).save()

    # lab A level
    lab_a_type = ln.Record(name="LabA", is_type=True).save()
    sample_a_type = ln.Record(name="Sample", is_type=True, type=lab_a_type).save()
    for name in ["s3", "s4"]:
        ln.Record(name=name, type=sample_a_type).save()

    # lab B level
    lab_b_type = ln.Record(name="LabB", is_type=True).save()
    sample_b_type = ln.Record(name="Sample", is_type=True, type=lab_b_type).save()
    for name in ["s5", "s6"]:
        ln.Record(name=name, type=sample_b_type).save()

    df = pd.DataFrame(
        {
            "biosample_name": pd.Categorical(["s1", "s2", "s3", "s4", "s5", "s6"]),
        }
    )

    feature = ln.Feature(name="biosample_name", dtype=sample_a_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert (
        "4 terms not validated in feature 'biosample_name': 's1', 's2', 's5', 's6'"
        in error.exconly()
    )
    assert set(curator.cat._cat_vectors["biosample_name"]._validated) == {
        "s3",
        "s4",
    }
    assert set(curator.cat._cat_vectors["biosample_name"]._non_validated) == {
        "s1",
        "s2",
        "s5",
        "s6",
    }

    # Move LabB under LabA
    lab_b_type.type = lab_a_type
    lab_b_type.save()
    feature.delete(permanent=True)  # re-create the feature with the new dtype
    feature = ln.Feature(name="biosample_name", dtype=lab_a_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert set(curator.cat._cat_vectors["biosample_name"]._validated) == {
        "s3",
        "s4",
        "s5",
        "s6",
    }
    assert set(curator.cat._cat_vectors["biosample_name"]._non_validated) == {
        "s1",
        "s2",
    }

    # Lab at the root
    feature.delete(permanent=True)  # re-create the feature with the new dtype
    feature = ln.Feature(name="biosample_name", dtype=sample_root_type).save()
    curator = ln.curators.DataFrameCurator(df, ln.examples.schemas.valid_features())
    with pytest.raises(ln.errors.ValidationError) as error:
        curator.validate()
    assert set(curator.cat._cat_vectors["biosample_name"]._validated) == {
        "s1",
        "s2",
    }
    assert set(curator.cat._cat_vectors["biosample_name"]._non_validated) == {
        "s3",
        "s4",
        "s5",
        "s6",
    }

    sample_a_type.records.all().delete(permanent=True)
    sample_b_type.records.all().delete(permanent=True)
    lab_b_type.records.all().delete(permanent=True)
    lab_a_type.records.all().delete(permanent=True)
    lab_a_type.delete(permanet=True)
    lab_b_type.delete(permanent=True)
    sample_root_type.records.all().delete(permanent=True)
    sample_root_type.delete(permanent=True)
    feature.delete(permanent=True)
