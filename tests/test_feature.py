import pandas as pd
from pandas.api.types import is_categorical_dtype, is_string_dtype

import lamindb as ln

df = pd.DataFrame(
    {
        "feat1": [1, 2, 3],
        "feat2": [3.1, 4.2, 5.3],
        "feat3": ["cond1", "cond2", "cond2"],
        "feat4": ["id1", "id2", "id3"],
    }
)


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
            assert feature.type == "Category"
        else:
            orig_type = df[feature.name].dtype.name
            orig_type_stripped = "".join(i for i in orig_type if not i.isdigit())
            assert feature.type == orig_type_stripped
    for feature in features:
        feature.save()
    assert set(ln.Category.select(feature__name="feat3").list("name")) == set(
        ["cond1", "cond2"]
    )
    for name in df.columns:
        queried_feature = ln.Feature.select(name=name).one()
        if name in categoricals:
            assert queried_feature.type == "Category"
        else:
            orig_type = df[name].dtype.name
            orig_type_stripped = "".join(i for i in orig_type if not i.isdigit())
            assert queried_feature.type == orig_type_stripped
    assert set(
        ln.Feature.select(name="feat3")
        .one()
        .categories.all()
        .values_list("name", flat=True)
    ) == set(["cond1", "cond2"])
    for feature in features:
        feature.delete()
    assert len(ln.Category.select(feature__name="feat3").list("name")) == 0
