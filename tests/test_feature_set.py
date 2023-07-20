from inspect import signature

import pandas as pd
import pytest

import lamindb as ln
from lamindb import _feature_set

df = pd.DataFrame(
    {
        "feat1": [1, 2, 3],
        "feat2": [3.1, 4.2, 5.3],
        "feat3": ["cond1", "cond2", "cond2"],
        "feat4": ["id1", "id2", "id3"],
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
    class_methods = ["from_values", "from_df"]
    for name in class_methods:
        setattr(Mock, name, getattr(_feature_set, name))
        assert signature(getattr(Mock, name)) == _feature_set.SIGS.pop(name)
    # methods
    for name, sig in _feature_set.SIGS.items():
        assert signature(getattr(_feature_set, name)) == sig


def test_feature_set_from_values():
    features = ["feat1", "feat2"]
    feature_set = ln.FeatureSet.from_values(features)
    id = feature_set.id
    assert feature_set._state.adding
    assert feature_set.type == "Feature"
    assert feature_set.schema == "core"
    assert feature_set.field == "name"
    feature_set.save()
    # test that the feature_set is retrieved from the database
    # in case it already exists
    feature_set = ln.FeatureSet.from_values(features)
    assert not feature_set._state.adding
    assert id == feature_set.id
    feature_set.delete()

    # edge cases
    with pytest.raises(ValueError):
        feature_set = ln.FeatureSet.from_values([])


def test_feature_set_from_records():
    features = ln.Feature.from_values(["feat1", "feat2"], "name")
    feature_set = ln.FeatureSet(features)
    id = feature_set.id
    assert feature_set._state.adding
    assert feature_set.type == "Feature"
    assert feature_set.schema == "core"
    assert feature_set.field == "id"
    feature_set.save()
    # test that the feature_set is retrieved from the database
    # in case it already exists
    feature_set = ln.FeatureSet(features)
    assert not feature_set._state.adding
    assert id == feature_set.id
    feature_set.delete()


def test_feature_set_from_df():
    feature_set = ln.FeatureSet.from_df(df)
    feature_set.save()
    assert set(ln.Category.select(feature__name="feat3").list("name")) == set(
        ["cond1", "cond2"]
    )
    for feature in feature_set.features.all():
        feature.delete()
    assert len(ln.Category.select(feature__name="feat3").list("name")) == 0
    feature_set.delete()
