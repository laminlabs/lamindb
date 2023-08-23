from inspect import signature

import lnschema_bionty as lb
import pandas as pd
import pytest
from lnschema_core.models import FileLabel
from pandas.api.types import is_categorical_dtype, is_string_dtype

import lamindb as ln
from lamindb import _feature

df = pd.DataFrame(
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


def test_feature_from_df():
    # try to generate the file without validated features
    feat1 = ln.Feature.filter(name="feat1").one_or_none()
    if feat1 is not None:
        feat1.delete()
    file = ln.File.from_df(df, description="test")
    assert file._feature_sets == {}
    # now, register all 4 features
    ln.save(ln.Feature.from_df(df.iloc[:, :4]))
    # try again
    file = ln.File.from_df(df, description="test")
    assert "columns" in file._feature_sets
    file.save()
    feature_set = file._feature_sets["columns"]
    features = feature_set.features.all()
    assert len(features) == len(df.columns[:4])
    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {col: df[col] for col in df.columns if is_categorical_dtype(df[col])}
    for key in string_cols:
        c = pd.Categorical(df[key])
        if len(c.categories) < len(c):
            categoricals[key] = c
    for feature in features:
        if feature.name in categoricals:
            assert feature.type == "category"
        else:
            orig_type = df[feature.name].dtype.name
            orig_type_stripped = "".join(i for i in orig_type if not i.isdigit())
            assert feature.type == orig_type_stripped
    for feature in features:
        feature.save()
    labels = [ln.Label(name=name) for name in df["feat3"].unique()]
    ln.save(labels)
    file.add_labels(labels, feature="feat3")
    assert set(ln.Label.filter(filelabel__feature__name="feat3").list("name")) == set(
        ["cond1", "cond2"]
    )
    for name in df.columns[:4]:
        queried_feature = ln.Feature.filter(name=name).one()
        if name in categoricals:
            assert queried_feature.type == "category"
        else:
            orig_type = df[name].dtype.name
            orig_type_stripped = "".join(i for i in orig_type if not i.isdigit())
            assert queried_feature.type == orig_type_stripped
    filelabel_links = FileLabel.objects.filter(file_id=file.id, feature__name="feat3")
    label_ids = filelabel_links.values_list("label_id")
    assert set(
        ln.Label.objects.filter(id__in=label_ids).values_list("name", flat=True)
    ) == set(["cond1", "cond2"])
    for feature in features:
        feature.delete()
    file.delete(storage=True)


def test_feature_init():
    # no args allowed
    with pytest.raises(ValueError):
        ln.Feature("x")
    # registries has to be a list of Registry types
    with pytest.raises(ValueError):
        ln.Feature(name="feat", type="category", registries=1)
    # each element of the list has to be a Registry
    with pytest.raises(ValueError):
        ln.Feature(name="feat", type="category", registries=[1])
    # registries_str
    feature = ln.Feature(name="feat1", type="category", registries=[ln.Label, lb.Gene])
    assert feature.registries == "core.Label|bionty.Gene"
