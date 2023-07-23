from inspect import signature

import pandas as pd
from pandas.api.types import is_categorical_dtype, is_string_dtype

import lamindb as ln
from lamindb import _feature

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
    class_methods = ["from_df"]
    for name in class_methods:
        setattr(Mock, name, getattr(_feature, name))
        assert signature(getattr(Mock, name)) == _feature.SIGS.pop(name)
    # methods
    for name, sig in _feature.SIGS.items():
        assert signature(getattr(_feature, name)) == sig


def test_feature_from_df():
    features = ln.Feature.from_df(df)
    assert len(features) == len(df.columns)
    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {col: df[col] for col in df.columns if is_categorical_dtype(df[col])}
    for key in string_cols:
        c = pd.Categorical(df[key])
        if len(c.categories) < len(c):
            categoricals[key] = c
    for feature in features:
        if feature.name in categoricals:
            assert feature.type == "categorical"
        else:
            orig_type = df[feature.name].dtype.name
            orig_type_stripped = "".join(i for i in orig_type if not i.isdigit())
            assert feature.type == orig_type_stripped
    for feature in features:
        feature.save()
    labels = ln.Label.from_values(df["feat3"], feature="feat3")
    ln.save(labels)
    assert set(ln.Label.select(feature__name="feat3").list("name")) == set(
        ["cond1", "cond2"]
    )
    for name in df.columns:
        queried_feature = ln.Feature.select(name=name).one()
        if name in categoricals:
            assert queried_feature.type == "categorical"
        else:
            orig_type = df[name].dtype.name
            orig_type_stripped = "".join(i for i in orig_type if not i.isdigit())
            assert queried_feature.type == orig_type_stripped
    assert set(
        ln.Feature.select(name="feat3")
        .one()
        .labels.all()
        .values_list("name", flat=True)
    ) == set(["cond1", "cond2"])
    for feature in features:
        feature.delete()
    assert len(ln.Label.select(feature__name="feat3").list("name")) == 0
