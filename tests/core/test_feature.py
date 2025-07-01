import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.errors import ValidationError
from lamindb.models.feature import serialize_pandas_dtype
from pandas.api.types import is_string_dtype


@pytest.fixture(scope="module")
def df():
    return pd.DataFrame(
        {
            "feat1": [1, 2, 3],
            "feat2": [3.1, 4.2, 5.3],
            "feat3": ["cond1", "cond2", "cond2"],
            "feat4": ["id1", "id2", "id3"],
            "rando_feature": ["rando1", "rando2", "rando3"],
        }
    )


def test_feature_init():
    # positional args not supported
    with pytest.raises(ValueError):
        ln.Feature("x")

    # dtype required unless is_type=True
    with pytest.raises(ValidationError):
        ln.Feature(name="feat")

    # is OK if also is_type is passed
    ln.Feature(name="Feat", is_type=True)

    # invalid dtype string
    with pytest.raises(ValueError):
        ln.Feature(name="feat", dtype="x")

    # categorical dtype must specify valid types
    with pytest.raises(ValidationError):
        ln.Feature(name="feat", dtype="cat[1]")

    # ensure feat1 does not exist
    if feat1 := ln.Feature.filter(name="feat1").one_or_none() is not None:
        feat1.delete()

    feat1 = ln.Feature(name="feat", dtype="str").save()
    # duplicate name with different dtype should fail
    with pytest.raises(ValidationError) as error:
        ln.Feature(name="feat", dtype="cat")
    assert (
        error.exconly()
        == "lamindb.errors.ValidationError: Feature feat already exists with dtype str, you passed cat"
    )
    feat1.delete()

    # string and list syntax for categorical dtypes should be equivalent and work
    feat2 = ln.Feature(name="feat2", dtype="str", description="feat2").save()
    feat2_again = ln.Feature(name="feat2", dtype="str", description="feat2").save()
    assert feat2 == feat2_again
    feat2.delete()

    # check that this works
    feature = ln.Feature(name="feat1", dtype="cat[ULabel|bionty.Gene]")
    # check that it also works via objects
    feature = ln.Feature(name="feat1", dtype=[ln.ULabel, bt.Gene])
    assert feature.dtype == "cat[ULabel|bionty.Gene]"


def test_cat_filters_dtype():
    feature = ln.Feature(
        name="disease",
        dtype=bt.Disease,
        cat_filters={
            "source__uid": "4a3ejKuf"
        },  # uid corresponds to disease_ontology_old.uid
    ).save()

    assert feature.dtype == "cat[bionty.Disease[source__uid='4a3ejKuf']]"

    feature.delete()


def test_cat_filters_empty_filter():
    # empty filter values should be rejected
    with pytest.raises(ValidationError) as error:
        ln.Feature(name="feat_empty", dtype=bt.Disease, cat_filters={"source__uid": ""})
    assert (
        "lamindb.errors.ValidationError: Empty value in filter source__uid"
        in error.exconly()
    )


def test_cat_filters_invalid_field_name():
    # invalid filter field names should be rejected
    source = bt.Source(
        name="", description="", organism="", entity="", version=""
    ).save()
    with pytest.raises(ValidationError) as error:
        ln.Feature(
            name="feat_invalid_attr",
            dtype=bt.Disease,
            cat_filters={"source__invalid_field": source},
        )
    assert (
        "lamindb.errors.ValidationError: SQLRecord Source has no attribute 'invalid_field' in filter source__invalid_field"
        in error.exconly()
    )
    source.delete()


def test_feature_from_df(df):
    if feat1 := ln.Feature.filter(name="feat1").one_or_none() is not None:
        feat1.delete()
    features = ln.Feature.from_df(df.iloc[:, :4]).save()
    artifact = ln.Artifact.from_df(df, description="test").save()
    # test for deprecated add_feature_set
    artifact.features._add_schema(ln.Schema(features), slot="columns")
    features = artifact.features.slots["columns"].features.all()
    assert len(features) == len(df.columns[:4])
    [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {
        col: df[col] for col in df.columns if isinstance(df[col], pd.CategoricalDtype)
    }
    for feature in features:
        if feature.name in categoricals:
            assert feature.dtype == "cat"
        else:
            orig_type = df[feature.name].dtype
            assert feature.dtype == serialize_pandas_dtype(orig_type)
    for feature in features:
        feature.save()
    labels = [ln.ULabel(name=name) for name in df["feat3"].unique()]
    ln.save(labels)
    feature = ln.Feature.get(name="feat3")
    feature.dtype = "cat"
    feature.save()
    with pytest.raises(ValidationError) as err:
        artifact.labels.add(labels, feature=feature)
    assert (
        err.exconly()
        == "lamindb.errors.ValidationError: Cannot manually annotate a feature measured *within* the dataset. Please use a Curator."
    )
    extfeature = ln.Feature(name="extfeat", dtype="str").save()
    with pytest.raises(ValidationError) as err:
        artifact.labels.add(labels, feature=extfeature)
    assert (
        err.exconly()
        == f"lamindb.errors.ValidationError: Feature {extfeature.name} needs dtype='cat' for label annotation, currently has dtype='str'"
    )

    # clean up
    artifact.delete(permanent=True)
    ln.Schema.filter().all().delete()
    ln.ULabel.filter().all().delete()
    ln.Feature.filter().all().delete()
