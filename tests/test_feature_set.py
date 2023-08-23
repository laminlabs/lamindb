from inspect import signature

import lnschema_bionty as lb
import pandas as pd
import pytest

import lamindb as ln
from lamindb import _feature_set
from lamindb._feature_set import get_related_name, sanity_check_features

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
    gene_symbols = ["TCF7", "MYC"]
    lb.settings.species = "human"
    feature_set = ln.FeatureSet.from_values(gene_symbols, lb.Gene.symbol)
    id = feature_set.id
    assert feature_set._state.adding
    assert feature_set.type == "float"
    assert feature_set.registry == "bionty.Gene"
    feature_set.save()
    # test that the feature_set is retrieved from the database
    # in case it already exists
    feature_set = ln.FeatureSet.from_values(gene_symbols, lb.Gene.symbol)
    assert not feature_set._state.adding
    assert id == feature_set.id
    feature_set.delete()

    # edge cases
    with pytest.raises(ValueError):
        feature_set = ln.FeatureSet.from_values([])

    with pytest.raises(TypeError):
        ln.FeatureSet.from_values(["a"], field="name")
    with pytest.raises(ValueError):
        ln.FeatureSet.from_values(["name"], field=ln.Feature.name, type="rna")
    with pytest.raises(TypeError):
        ln.FeatureSet.from_values([1], field=ln.Label.name, type="rna")

    # return none if no validated features
    assert ln.FeatureSet.from_values(["name"], field=ln.Label.name, type="rna") is None


def test_feature_set_from_records():
    features = ln.Feature.from_df(df)
    feature_set = ln.FeatureSet(features)
    id = feature_set.id
    assert feature_set._state.adding
    assert feature_set.type is None
    assert feature_set.registry == "core.Feature"
    feature_set.save()
    # test that the feature_set is retrieved from the database
    # in case it already exists
    feature_set = ln.FeatureSet(features)
    assert not feature_set._state.adding
    assert id == feature_set.id
    feature_set.delete()

    # edge case
    with pytest.raises(ValueError):
        positional_arg = 1
        ln.FeatureSet(features, positional_arg)


def test_feature_set_from_df():
    feature_set = ln.FeatureSet.from_df(df)
    feature_set.save()
    assert feature_set.type is None
    for feature in feature_set.features.all():
        feature.delete()
    feature_set.delete()


def test_get_related_name():
    with pytest.raises(ValueError):
        get_related_name(ln.Transform)


def test_sanity_check_features():
    with pytest.raises(ValueError):
        sanity_check_features([])
    with pytest.raises(TypeError):
        sanity_check_features(["feature"])
    with pytest.raises(TypeError):
        sanity_check_features({"feature"})
    with pytest.raises(ValueError):
        sanity_check_features([ln.Run(), ln.Transform()])


def test_kwargs():
    with pytest.raises(ValueError):
        ln.FeatureSet(x="1", features=[])


def test_pass_modality():
    m = ln.Modality(name="rna")
    feature_set = ln.FeatureSet(modality=m, features=[ln.Feature(name="f", type="f")])
    assert feature_set.modality == m
    with pytest.raises(ValueError):
        feature_set = ln.FeatureSet(
            modality=ln.Transform(), features=[ln.Feature(name="f", type="f")]
        )
