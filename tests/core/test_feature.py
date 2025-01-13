from inspect import signature

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb import _feature
from lamindb._feature import convert_pandas_dtype_to_lamin_dtype
from lamindb.core.exceptions import ValidationError
from pandas.api.types import is_categorical_dtype, is_string_dtype


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


def test_signatures():
    # this seems currently the easiest and most transparent
    # way to test violations of the signature equality
    # the MockORM class is needed to get inspect.signature
    # to work
    class Mock:
        pass

    # class methods
    class_methods = ["from_df"]
    for name in class_methods:
        setattr(Mock, name, getattr(_feature, name))
        assert signature(getattr(Mock, name)) == _feature.SIGS.pop(name)
    # methods
    for name, sig in _feature.SIGS.items():
        assert signature(getattr(_feature, name)) == sig


def test_feature_init():
    # no args allowed
    with pytest.raises(ValueError):
        ln.Feature("x")
    # no type passed
    with pytest.raises(ValueError):
        ln.Feature(name="feat")
    # wrong type
    with pytest.raises(ValueError):
        ln.Feature(name="feat", dtype="x")
    # type has to be a list of Record types
    with pytest.raises(ValueError):
        ln.Feature(name="feat", dtype="cat[1]")
    # ensure feat1 does not exist
    if feat1 := ln.Feature.filter(name="feat1").one_or_none() is not None:
        feat1.delete()
    feat1 = ln.Feature(name="feat", dtype="str").save()
    with pytest.raises(ValidationError) as error:
        ln.Feature(name="feat", dtype="cat")
    assert (
        error.exconly()
        == "lamindb.core.exceptions.ValidationError: Feature feat already exists with dtype str, you passed cat"
    )
    feat1.delete()
    # check that this works
    feature = ln.Feature(name="feat1", dtype="cat[ULabel|bionty.Gene]")
    # check that it also works via objects
    feature = ln.Feature(name="feat1", dtype=[ln.ULabel, bt.Gene])
    assert feature.dtype == "cat[ULabel|bionty.Gene]"


def test_feature_from_df(df):
    if feat1 := ln.Feature.filter(name="feat1").one_or_none() is not None:
        feat1.delete()
    features = ln.Feature.from_df(df.iloc[:, :4]).save()
    artifact = ln.Artifact.from_df(df, description="test").save()
    artifact.features.add_schema(ln.Schema(features), slot="columns")
    features = artifact.features["columns"]
    assert len(features) == len(df.columns[:4])
    [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {col: df[col] for col in df.columns if is_categorical_dtype(df[col])}
    for feature in features:
        if feature.name in categoricals:
            assert feature.dtype == "cat"
        else:
            orig_type = df[feature.name].dtype
            assert feature.dtype == convert_pandas_dtype_to_lamin_dtype(orig_type)
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
        == "lamindb.core.exceptions.ValidationError: Cannot manually annotate internal feature with label. Please use ln.Curator"
    )
    extfeature = ln.Feature(name="extfeat", dtype="str").save()
    with pytest.raises(ValidationError) as err:
        artifact.labels.add(labels, feature=extfeature)
    assert (
        err.exconly()
        == f"lamindb.core.exceptions.ValidationError: Feature {extfeature.name} needs dtype='cat' for label annotation, currently has dtype='str'"
    )

    # clean up
    artifact.delete(permanent=True)
    ln.Schema.filter().all().delete()
    ln.ULabel.filter().all().delete()
    ln.Feature.filter().all().delete()
