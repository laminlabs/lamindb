from inspect import signature

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb import _schema
from lamindb._schema import get_related_name, validate_features
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
        setattr(Mock, name, getattr(_schema, name))
        assert signature(getattr(Mock, name)) == _schema.SIGS.pop(name)
    # methods
    for name, sig in _schema.SIGS.items():
        assert signature(getattr(_schema, name)) == sig


def test_schema_from_values():
    gene_symbols = ["TCF7", "MYC"]
    bt.settings.organism = "human"
    bt.Gene.filter(symbol__in=gene_symbols).all().delete()
    with pytest.raises(ValidationError) as error:
        schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert error.exconly().startswith(
        "lamindb.core.exceptions.ValidationError: These values could not be validated:"
    )
    ln.save(bt.Gene.from_values(gene_symbols, "symbol"))
    schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol)
    # below should be a queryset and not a list
    assert set(schema.members) == set(bt.Gene.from_values(gene_symbols, "symbol"))
    assert schema.dtype == "num"  # this is NUMBER_TYPE
    schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert schema._state.adding
    assert schema.dtype == "int"
    assert schema.registry == "bionty.Gene"
    schema.save()
    assert set(schema.members) == set(schema.genes.all())
    id = schema.id
    # test that the schema is retrieved from the database
    # in case it already exists
    schema = ln.Schema.from_values(gene_symbols, bt.Gene.symbol, type=int)
    assert not schema._state.adding
    assert id == schema.id
    schema.delete()

    # edge cases
    with pytest.raises(ValueError):
        schema = ln.Schema.from_values([])

    with pytest.raises(TypeError):
        ln.Schema.from_values(["a"], field="name")
    with pytest.raises(ValidationError):
        schema = ln.Schema.from_values(
            ["weird_name"], field=ln.Feature.name, type="float"
        )
    with pytest.raises(ValidationError):
        ln.Schema.from_values([1], field=ln.Feature.name, type="float")

    # return none if no validated features
    with pytest.raises(ValidationError):
        ln.Schema.from_values(["name"], field=ln.Feature.name, type="float")


def test_schema_from_records(df):
    features = ln.Feature.from_df(df)
    with pytest.raises(ValueError) as error:
        schema = ln.Schema(features)
    assert (
        error.exconly()
        == "ValueError: Can only construct feature sets from validated features"
    )

    ln.save(features)
    schema = ln.Schema(features)
    assert schema.id is None
    assert schema._state.adding
    assert schema.dtype is None
    assert schema.registry == "Feature"
    schema.save()
    # test that the schema is retrieved from the database
    # in case it already exists
    schema = ln.Schema(features)
    assert not schema._state.adding
    assert schema.id is not None
    schema.delete()

    # edge case
    with pytest.raises(ValueError):
        positional_arg = 1
        ln.Schema(features, positional_arg)


def test_schema_from_df(df):
    # test using type
    bt.settings.organism = "human"
    genes = [bt.Gene(symbol=name) for name in df.columns]
    ln.save(genes)
    with pytest.raises(ValueError) as error:
        ln.Schema.from_df(df, field=bt.Gene.symbol)
    assert error.exconly().startswith("ValueError: data types are heterogeneous:")
    schema = ln.Schema.from_df(df[["feat1", "feat2"]], field=bt.Gene.symbol)
    for gene in genes:
        gene.delete()

    # now for the features registry
    features = ln.Feature.from_df(df)
    ln.save(features)
    schema = ln.Schema.from_df(df).save()
    assert schema.dtype is None
    ln.Schema.filter().all().delete()
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
    transform = ln.Transform(key="test")
    transform.save()
    # This is just a type check
    with pytest.raises(TypeError) as error:
        validate_features([transform, ln.Run(transform)])
    assert error.exconly() == "TypeError: schema can only contain a single type"
    transform.delete()


def test_kwargs():
    with pytest.raises(ValueError):
        ln.Schema(x="1", features=[])


def test_edge_cases():
    feature = ln.Feature(name="rna", dtype="float")
    ln.save([feature])
    with pytest.raises(ValueError) as error:
        ln.Schema(feature)
    assert (
        error.exconly()
        == "ValueError: Please pass a ListLike of features, not a single feature"
    )
    feature.delete()
