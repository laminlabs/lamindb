from inspect import signature

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb import _feature
from lamindb._feature import convert_numpy_dtype_to_lamin_feature_type
from lnschema_core.models import ArtifactULabel
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


def test_feature_from_df(df):
    # try to generate the file without validated features
    feat1 = ln.Feature.filter(name="feat1").one_or_none()
    if feat1 is not None:
        feat1.delete()
    artifact = ln.Artifact.from_df(df, description="test")

    # now, register all 4 features
    features = ln.Feature.from_df(df.iloc[:, :4])
    ln.save(features)
    # try again
    artifact = ln.Artifact.from_df(df, description="test")
    artifact.save()
    # link features
    artifact.features.add_feature_set(ln.FeatureSet(features), slot="columns")
    features = artifact.features["columns"]
    assert len(features) == len(df.columns[:4])
    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {col: df[col] for col in df.columns if is_categorical_dtype(df[col])}
    for key in string_cols:
        c = pd.Categorical(df[key])
        if len(c.categories) < len(c):
            categoricals[key] = c
    for feature in features:
        if feature.name in categoricals:
            assert feature.dtype == "cat"
        else:
            orig_type = df[feature.name].dtype
            assert feature.dtype == convert_numpy_dtype_to_lamin_feature_type(orig_type)
    for feature in features:
        feature.save()
    labels = [ln.ULabel(name=name) for name in df["feat3"].unique()]
    ln.save(labels)
    features_lookup = ln.Feature.lookup()
    artifact.labels.add(labels, feature=features_lookup.feat3)
    assert set(
        ln.ULabel.filter(artifact_links__feature__name="feat3").list("name")
    ) == {"cond1", "cond2"}
    for name in df.columns[:4]:
        queried_feature = ln.Feature.filter(name=name).one()
        if name in categoricals:
            assert queried_feature.dtype == "cat[ULabel]"
        else:
            orig_type = df[name].dtype
            assert queried_feature.dtype == convert_numpy_dtype_to_lamin_feature_type(
                orig_type
            )
    links_artifactlabel = ArtifactULabel.objects.filter(
        artifact_id=artifact.id, feature__name="feat3"
    )
    label_ids = links_artifactlabel.values_list("ulabel_id")
    assert set(
        ln.ULabel.objects.filter(id__in=label_ids).values_list("name", flat=True)
    ) == {"cond1", "cond2"}

    # clean up
    artifact.delete(permanent=True)
    ln.FeatureSet.filter().all().delete()
    ln.ULabel.filter().all().delete()
    ln.Feature.filter().all().delete()


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
    feat1 = ln.Feature.filter(name="feat1").one_or_none()
    if feat1 is not None:
        feat1.delete()
    # check that this works
    feature = ln.Feature(name="feat1", dtype="cat[ULabel|bionty.Gene]")
    # check that it also works via objects
    feature = ln.Feature(name="feat1", dtype=[ln.ULabel, bt.Gene])
    assert feature.dtype == "cat[ULabel|bionty.Gene]"
