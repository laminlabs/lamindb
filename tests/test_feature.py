import pandas as pd
from pandas.api.types import is_categorical_dtype, is_string_dtype

import lamindb as ln

df = pd.DataFrame({"feat1": [1, 2], "feat2": [3.1, 4.2], "feat3": ["cond1", "cond2"]})


def test_feature_from_df():
    features = ln.Feature.from_df(df)
    assert len(features) == len(df.columns)
    string_or_categorical_columns = [
        col
        for col in df.columns
        if is_string_dtype(df[col]) or is_categorical_dtype(df[col])
    ]
    for feature in features:
        if feature.name in string_or_categorical_columns:
            assert feature.type == "str"
        else:
            assert feature.type == df[feature.name].dtype.name
    for feature in features:
        feature.save()
    assert set(ln.FeatureValue.select(feature__name="feat3").list("value")) == set(
        ["cond1", "cond2"]
    )
    for name in df.columns:
        queried_feature = ln.Feature.select(name=name).one()
        if name in string_or_categorical_columns:
            assert queried_feature.type == "str"
        else:
            assert queried_feature.type == df[name].dtype.name
    assert set(
        ln.Feature.select(name="feat3")
        .one()
        .values.all()
        .values_list("value", flat=True)
    ) == set(["cond1", "cond2"])
    for feature in features:
        feature.delete()
    assert len(ln.FeatureValue.select(feature__name="feat3").list("value")) == 0
