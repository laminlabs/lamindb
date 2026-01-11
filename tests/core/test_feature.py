import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamindb.errors import ValidationError
from lamindb.models.feature import serialize_pandas_dtype
from pandas.api.types import is_string_dtype


@pytest.fixture(scope="module")
def dict_data():
    return {
        "dict_feat1": 42,
        "dict_feat2": 3.14,
        "dict_feat3": "somestring",  # string (ambiguous cat ? str)
        "dict_feat4": True,
        "dict_feat5": [1, 2, 3],
        "dict_feat6": ["a", "b", "c"],  # list[str] (ambiguous list[cat ? str])
        "dict_feat7": {"key": "value"},
    }


def test_feature_init():
    # positional args not supported
    with pytest.raises(ValueError):
        ln.Feature("x")

    # dtype required unless is_type=True
    with pytest.raises(ValidationError):
        ln.Feature(name="feat")

    # is OK if also is_type is passed
    ln.Feature(name="Feat", is_type=True)

    # invalid dtype string
    with pytest.raises(ValueError):
        ln.Feature(name="feat", dtype="x")

    # categorical dtype must specify valid types
    with pytest.raises(ValidationError):
        ln.Feature(name="feat", dtype="cat[1]")

    # ensure feat1 does not exist
    if feat1 := ln.Feature.filter(name="feat1").one_or_none() is not None:
        feat1.delete(permanent=True)

    feat1 = ln.Feature(name="feat", dtype="str").save()
    # duplicate name with different dtype should fail
    with pytest.raises(ValidationError) as error:
        ln.Feature(name="feat", dtype=ln.ULabel)
    assert (
        error.exconly()
        == "lamindb.errors.ValidationError: Feature feat already exists with dtype str, you passed cat[ULabel]"
    )
    feat1.delete(permanent=True)

    # string and list syntax for categorical dtypes should be equivalent and work
    feat2 = ln.Feature(name="feat2", dtype="str", description="feat2").save()
    feat2_again = ln.Feature(name="feat2", dtype="str", description="feat2").save()
    assert feat2 == feat2_again
    feat2.delete(permanent=True)

    # categorical dtype with union of registries using string syntax must be valid
    feature = ln.Feature(name="feat1", dtype="cat[Record|bionty.Gene]")
    assert feature._dtype_str == "cat[Record|bionty.Gene]"
    # categorical dtype with union of registries using objects must be valid
    feature = ln.Feature(name="feat1", dtype=[ln.Record, bt.Gene])
    assert feature._dtype_str == "cat[Record|bionty.Gene]"

    # dtype with field name before bracket filters must be valid
    feature = ln.Feature(
        name="gene_feature", dtype="cat[bionty.Gene.ensembl_gene_id[organism='human']]"
    )
    print(feature._dtype_str)
    assert "bionty.Gene" in feature._dtype_str
    assert "ensembl_gene_id" in feature._dtype_str
    assert "organism='human'" in feature._dtype_str


# @pytest.mark.skipif(
#     os.getenv("LAMINDB_TEST_DB_VENDOR") == "sqlite", reason="Postgres-only"
# )
# def test_cannot_mutate_dtype():
#     feature = ln.Feature(name="feature", dtype=str).save()
#     feature._dtype_str = int
#     with pytest.raises(django.db.utils.IntegrityError) as error:
#         feature.save()
#     assert "dtype field is immutable and cannot be changed" in error.exconly()
#     feature.delete(permanent=True)


# def test_cat_filters_dtype():
#     feature = ln.Feature(
#         name="disease",
#         dtype=bt.Disease,
#         cat_filters={
#             "source__uid": "4a3ejKuf"
#         },  # uid corresponds to disease_ontology_old.uid
#     ).save()

#     assert feature._dtype_str == "cat[bionty.Disease[source__uid='4a3ejKuf']]"

#     feature.delete(permanent=True)


def test_cat_filters_empty_filter():
    # empty filter values should be rejected
    with pytest.raises(ValidationError) as error:
        ln.Feature(name="feat_empty", dtype=bt.Disease, cat_filters={"source__uid": ""})
    assert (
        "lamindb.errors.ValidationError: Empty value in filter source__uid"
        in error.exconly()
    )


def test_cat_filters_invalid_field_name():
    # invalid filter field names should be rejected
    source = bt.Source(
        name="", description="", organism="", entity="", version=""
    ).save()
    with pytest.raises(ValidationError) as error:
        ln.Feature(
            name="feat_invalid_attr",
            dtype=bt.Disease,
            cat_filters={"source__invalid_field": source},
        )
    assert (
        "lamindb.errors.ValidationError: SQLRecord Source has no attribute 'invalid_field' in filter source__invalid_field"
        in error.exconly()
    )
    source.delete(permanent=True)


def test_feature_from_df():
    df = pd.DataFrame(
        {
            "feat1": [1, 2, 3],
            "feat2": [3.1, 4.2, 5.3],
            "feat3": pd.Categorical(["cond1", "cond2", "cond2"]),
            "feat4": ["id1", "id2", "id3"],
            "rando_feature": ["rando1", "rando2", "rando3"],
        }
    )
    if feat1 := ln.Feature.filter(name="feat1").one_or_none() is not None:
        feat1.delete(permanent=True)
    features = ln.Feature.from_dataframe(df.iloc[:, :4]).save()
    artifact = ln.Artifact.from_dataframe(df, description="test").save()
    # test for deprecated add_feature_set
    schema = ln.Schema(features).save()
    artifact.features._add_schema(schema, slot="columns")
    features = artifact.features.slots["columns"].features.all()
    assert len(features) == len(df.columns[:4])
    [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {
        col: df[col] for col in df.columns if isinstance(df[col], pd.CategoricalDtype)
    }
    for feature in features:
        if feature.name in categoricals:
            assert feature._dtype_str == "cat"
        else:
            orig_type = df[feature.name].dtype
            assert feature._dtype_str == serialize_pandas_dtype(orig_type)
    for feature in features:
        feature.save()
    labels = [ln.Record(name=name) for name in df["feat3"].unique()]
    ln.save(labels)
    feature = ln.Feature.get(name="feat3")
    with pytest.raises(ValidationError) as err:
        artifact.labels.add(labels, feature=feature)
    assert (
        err.exconly()
        == "lamindb.errors.ValidationError: Cannot manually annotate a feature measured *within* the dataset. Please use a Curator."
    )
    extfeature = ln.Feature(name="extfeat", dtype="str").save()
    with pytest.raises(ValidationError) as err:
        artifact.labels.add(labels, feature=extfeature)
    assert (
        err.exconly()
        == f"lamindb.errors.ValidationError: Feature {extfeature.name} needs dtype='cat' for label annotation, currently has dtype='str'"
    )

    # clean up
    artifact.delete(permanent=True)
    ln.Schema.filter().delete(permanent=True)
    ln.Record.filter().delete(permanent=True)
    ln.Feature.filter().delete(permanent=True)


def test_feature_from_dict(dict_data):
    # defaults to str for ambiguous types
    features = ln.Feature.from_dict(dict_data)
    assert len(features) == len(dict_data)
    assert features[0]._dtype_str == "int"
    assert features[1]._dtype_str == "float"
    assert features[2]._dtype_str == "str"
    assert features[3]._dtype_str == "bool"
    assert features[4]._dtype_str == "list[int]"
    assert features[5]._dtype_str == "list[str]"
    assert features[6]._dtype_str == "dict"

    # Wrong field
    with pytest.raises(ValueError) as e:
        ln.Feature.from_dict(dict_data, field=ln.Record.name)
    assert "field must be a Feature FieldAttr" in str(e.value)

    # Explicit field
    features_with_field = ln.Feature.from_dict(dict_data, field=ln.Feature.name)
    assert len(features_with_field) == len(dict_data)


def test_feature_from_dict_type(dict_data):
    feature_type = ln.Feature(name="Testdata_feature_type", is_type=True).save()
    features = ln.Feature.from_dict(dict_data, type=feature_type).save()
    for feature in features:
        assert feature.type.name == "Testdata_feature_type"
    ln.Feature.filter(type__isnull=False).delete(permanent=True)
    feature_type.delete(permanent=True)


def test_feature_query_by_dtype():
    """Test querying Feature by dtype (deprecated) and _dtype_str."""
    str_feat = ln.Feature(name="test_str_feat", dtype=str).save()
    int_feat = ln.Feature(name="test_int_feat", dtype=int).save()
    try:
        # Test querying by _dtype_str (current way)
        str_features = ln.Feature.filter(_dtype_str="str", name="test_str_feat")
        assert str_features.count() == 1
        assert str_features.first() == str_feat

        str_features = ln.Feature.filter(dtype_as_str="str", name="test_str_feat")
        assert str_features.count() == 1
        assert str_features.first() == str_feat

        # Test querying by dtype (deprecated) - should work but issue warning
        with pytest.warns(
            DeprecationWarning,
            match="Querying Feature by `dtype` is deprecated.*Notice the new dtype encoding format",
        ):
            str_features_deprecated = ln.Feature.filter(
                dtype="str", name="test_str_feat"
            )
            assert str_features_deprecated.count() == 1
            assert str_features_deprecated.first() == str_feat
    finally:
        # Clean up
        str_feat.delete(permanent=True)
        int_feat.delete(permanent=True)
