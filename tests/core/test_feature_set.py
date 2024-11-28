from inspect import signature

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb import _feature_set
from lamindb._feature_set import get_related_name, validate_features
from lamindb.core.exceptions import ValidationError


@pytest.fixture(scope="module")
def df():
    return pd.DataFrame(
        {
            "feat1": [1, 2, 3],
            "feat2": [3, 4, 5],
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
    bt.settings.organism = "human"
    bt.Gene.filter(symbol__in=gene_symbols).all().delete()
    with pytest.raises(ValidationError) as error:
        feature_set = ln.FeatureSet.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert error.exconly().startswith(
        "lamindb.core.exceptions.ValidationError: These values could not be validated:"
    )
    ln.save(bt.Gene.from_values(gene_symbols, "symbol"))
    feature_set = ln.FeatureSet.from_values(gene_symbols, bt.Gene.symbol)
    # below should be a queryset and not a list
    assert set(feature_set.members) == set(bt.Gene.from_values(gene_symbols, "symbol"))
    assert feature_set.dtype == "num"  # this is NUMBER_TYPE
    feature_set = ln.FeatureSet.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert feature_set._state.adding
    assert feature_set.dtype == "int"
    assert feature_set.registry == "bionty.Gene"
    feature_set.save()
    assert set(feature_set.members) == set(feature_set.genes.all())
    id = feature_set.id
    # test that the feature_set is retrieved from the database
    # in case it already exists
    feature_set = ln.FeatureSet.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert not feature_set._state.adding
    assert id == feature_set.id
    feature_set.delete()

    # edge cases
    with pytest.raises(ValueError):
        feature_set = ln.FeatureSet.from_values([])

    with pytest.raises(TypeError):
        ln.FeatureSet.from_values(["a"], field="name")
    with pytest.raises(ValidationError):
        feature_set = ln.FeatureSet.from_values(
            ["weird_name"], field=ln.Feature.name, type="float"
        )
    with pytest.raises(ValidationError):
        ln.FeatureSet.from_values([1], field=ln.ULabel.name, type="float")

    # return none if no validated features
    with pytest.raises(ValidationError):
        ln.FeatureSet.from_values(["name"], field=ln.ULabel.name, type="float")


def test_feature_set_from_records(df):
    features = ln.Feature.from_df(df)
    with pytest.raises(ValueError) as error:
        feature_set = ln.FeatureSet(features)
    assert (
        error.exconly()
        == "ValueError: Can only construct feature sets from validated features"
    )

    ln.save(features)
    feature_set = ln.FeatureSet(features)
    assert feature_set.id is None
    assert feature_set._state.adding
    assert feature_set.dtype is None
    assert feature_set.registry == "Feature"
    feature_set.save()
    # test that the feature_set is retrieved from the database
    # in case it already exists
    feature_set = ln.FeatureSet(features)
    assert not feature_set._state.adding
    assert feature_set.id is not None
    feature_set.delete()

    # edge case
    with pytest.raises(ValueError):
        positional_arg = 1
        ln.FeatureSet(features, positional_arg)


def test_feature_set_from_df(df):
    # test using type
    bt.settings.organism = "human"
    genes = [bt.Gene(symbol=name) for name in df.columns]
    ln.save(genes)
    with pytest.raises(ValueError) as error:
        ln.FeatureSet.from_df(df, field=bt.Gene.symbol)
    assert error.exconly().startswith("ValueError: data types are heterogeneous:")
    feature_set = ln.FeatureSet.from_df(df[["feat1", "feat2"]], field=bt.Gene.symbol)
    for gene in genes:
        gene.delete()

    # now for the features registry
    features = ln.Feature.from_df(df)
    ln.save(features)
    feature_set = ln.FeatureSet.from_df(df)
    feature_set.save()
    assert feature_set.dtype is None
    ln.FeatureSet.filter().all().delete()
    ln.Feature.filter().all().delete()


def test_get_related_name():
    with pytest.raises(ValueError):
        get_related_name(ln.Transform)


def test_validate_features():
    with pytest.raises(ValueError):
        validate_features([])
    with pytest.raises(TypeError):
        validate_features(["feature"])
    with pytest.raises(TypeError):
        validate_features({"feature"})
    transform = ln.Transform(name="test")
    transform.save()
    # This is just a type check
    with pytest.raises(TypeError) as error:
        validate_features([transform, ln.Run(transform)])
    assert error.exconly() == "TypeError: feature_set can only contain a single type"
    transform.delete()


def test_kwargs():
    with pytest.raises(ValueError):
        ln.FeatureSet(x="1", features=[])


def test_edge_cases():
    feature = ln.Feature(name="rna", dtype="float")
    ln.save([feature])
    with pytest.raises(ValueError) as error:
        ln.FeatureSet(feature)
    assert (
        error.exconly()
        == "ValueError: Please pass a ListLike of features, not a single feature"
    )
    feature.delete()
