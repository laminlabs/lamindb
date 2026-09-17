from datetime import datetime

import bionty as bt
import lamindb as ln
import pandas as pd
import pytest
from lamin_utils import logger
from lamindb.errors import FieldValidationError, ValidationError
from lamindb.models.feature import (
    FeaturePredicate,
    _format_cat_filter_value,
    _split_filter_parts,
    convert_to_pandas_dtype,
    dtype_as_object,
    serialize_dtype,
    serialize_pandas_dtype,
)
from lamindb.models.record import get_feature_sqlrecord_field
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
    type_feat = ln.Feature(name="Feat", is_type=True)
    with pytest.warns(
        DeprecationWarning,
        match="Use dtype_as_str instead of dtype",
    ):
        assert type_feat.dtype is None

    # invalid dtype string
    with pytest.raises(ValueError):
        ln.Feature(name="feat", dtype="x")

    # categorical dtype must specify valid types
    with pytest.raises(ValidationError):
        ln.Feature(name="feat", dtype="cat[1]")
    # deprecated `coerce_dtype` should warn and still set coerce
    with pytest.warns(
        DeprecationWarning,
        match="`coerce_dtype` argument was renamed to `coerce`",
    ):
        feature = ln.Feature(
            name="feat-coerce-deprecated",
            dtype="str",
            coerce_dtype=True,
        )
    assert feature.coerce is True
    with pytest.warns(
        DeprecationWarning,
        match="Use coerce instead of coerce_dtype",
    ):
        assert feature.coerce_dtype is True
    feature.coerce_dtype = False
    assert feature.coerce is False
    with pytest.warns(
        DeprecationWarning,
        match="Use dtype_as_str instead of dtype",
    ):
        assert feature.dtype == "str"
    # unknown keyword args should raise a field validation error
    with pytest.raises(FieldValidationError):
        ln.Feature(name="feat", dtype="str", not_a_valid_kwarg=True)

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

    # categorical dtype with union of registry fields using objects must be valid
    feature = ln.Feature(
        name="feat1", dtype=[bt.Tissue.ontology_id, bt.CellType.ontology_id]
    )
    assert (
        feature._dtype_str
        == "cat[bionty.Tissue.ontology_id|bionty.CellType.ontology_id]"
    )

    # dtype with field name before bracket filters must be valid
    feature = ln.Feature(
        name="gene_feature", dtype="cat[bionty.Gene.ensembl_gene_id[organism='human']]"
    )
    print(feature._dtype_str)
    assert "bionty.Gene" in feature._dtype_str
    assert "ensembl_gene_id" in feature._dtype_str
    assert "organism='human'" in feature._dtype_str


def test_feature_values_through_roundtrip():
    author_feature = ln.Feature(name="values-from-author", dtype=ln.Record)
    assert author_feature._aux is None
    assert author_feature._related_feature_uid is None
    author_feature.save()
    books_feature = ln.Feature(
        name="values-from-books",
        dtype=list[ln.Record],
        values_through=author_feature,
    ).save()
    try:
        assert books_feature._aux["vf"] == author_feature.uid
        assert books_feature.values_through.uid == author_feature.uid
        assert books_feature.related_feature.uid == author_feature.uid
        assert author_feature.related_feature.uid == books_feature.uid
        reloaded_books_feature = ln.Feature.get(uid=books_feature.uid)
        assert reloaded_books_feature.values_through.uid == author_feature.uid

        # Clearing values_through should remove both forward and reverse relation metadata.
        books_feature.values_through = None
        books_feature.save()
        books_feature.refresh_from_db()
        author_feature.refresh_from_db()
        assert books_feature.values_through is None
        assert books_feature.related_feature is None
        assert books_feature._aux is None or "vf" not in books_feature._aux
        assert author_feature.related_feature is None
        assert author_feature._aux is None or "rf" not in author_feature._aux
    finally:
        books_feature.delete(permanent=True)
        author_feature.delete(permanent=True)


def test_feature_values_through_sqlrecord_field_roundtrip():
    feature = ln.Feature(
        name="values-from-created-at",
        dtype="datetime64[ns, UTC]",
        values_through="created_at",
    ).save()
    try:
        assert feature._aux["sf"] == "created_at"
        assert feature._aux.get("vf") is None
        assert feature.values_through == "created_at"
        assert feature.related_feature is None
        reloaded = ln.Feature.get(uid=feature.uid)
        assert reloaded.values_through == "created_at"
        assert reloaded._aux["sf"] == "created_at"
        assert get_feature_sqlrecord_field(reloaded) == "created_at"

        feature.values_through = None
        feature.save()
        feature.refresh_from_db()
        assert feature.values_through is None
        assert feature._aux is None or "sf" not in feature._aux
        feature._aux = None
        assert feature._sqlrecord_field is None

        with pytest.raises(
            TypeError,
            match="Feature.values_through expects a Feature, SQLRecordFieldName, or None",
        ):
            feature.values_through = 1
        with pytest.raises(
            ValueError, match="Unsupported feature field mapping 'extra_data'"
        ):
            feature.values_through = "extra_data"

        feature.values_through = "created_at"
        feature.save()
        feature.refresh_from_db()
        assert feature.values_through == "created_at"
        assert feature._aux["sf"] == "created_at"
    finally:
        feature.delete(permanent=True)


def test_feature_values_through_requires_saved_source():
    unsaved_source = ln.Feature(name="values-from-unsaved-source", dtype=ln.Record)
    with pytest.raises(
        AssertionError,
        match="requires a saved Feature object",
    ):
        ln.Feature(
            name="values-from-unsaved-target",
            dtype=list[ln.Record],
            values_through=unsaved_source,
        ).save()


def test_feature_values_through_setter_requires_no_existing_links():
    target = ln.Feature(name="values-from-setter-target", dtype=list[ln.Record]).save()
    source = ln.Feature(name="values-from-setter-source", dtype=ln.Record).save()
    record_a = ln.Record(name="values-from-setter-a").save()
    record_b = ln.Record(name="values-from-setter-b").save()
    try:
        record_a.features.set_values({"values-from-setter-target": [record_b]})
        with pytest.raises(
            ValueError,
            match="can only be set when no RecordRecord links exist",
        ):
            target.values_through = source
    finally:
        record_a.delete(permanent=True)
        record_b.delete(permanent=True)
        target.delete(permanent=True)
        source.delete(permanent=True)


def test_feature_predicate_cannot_cast_to_bool():
    feature = ln.Feature(name="predicate-bool-guard", dtype="str")
    predicate = feature == "x"
    with pytest.raises(
        TypeError,
        match="Feature predicates cannot be used as booleans",
    ):
        bool(predicate)


def test_should_build_model_predicate_returns_false_for_type_features():
    feature_type = ln.Feature(name="predicate-type-feature", is_type=True)
    other_model = ln.Feature(name="predicate-other-model", dtype="str")
    assert feature_type._should_build_model_predicate(other_model) is False


def test_feature_predicate_model_ne_and_ordering_comparators():
    feature = ln.Feature(name="predicate-model-ne", dtype=ln.Record)
    record = ln.Record(name="predicate-model-ne-record")
    model_predicate = feature != record
    assert isinstance(model_predicate, FeaturePredicate)
    assert model_predicate.comparator == "__ne"
    assert model_predicate.value is record

    ge_predicate = feature >= 1
    lt_predicate = feature < 1
    assert isinstance(ge_predicate, FeaturePredicate)
    assert isinstance(lt_predicate, FeaturePredicate)
    assert ge_predicate.comparator == "__gte"
    assert lt_predicate.comparator == "__lt"


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


@pytest.mark.parametrize("filter_value", [None, "", [], 0])
def test_cat_filters_empty_filter(filter_value):
    # empty filter values should be rejected
    with pytest.raises(ValidationError) as error:
        ln.Feature(
            name="feat_empty",
            dtype=bt.Disease,
            cat_filters={"source__uid": filter_value},
        )
    assert "Empty value in filter source__uid" in error.exconly()


def test_cat_filters_incompatible_with_nested_dtype():
    with pytest.raises(ValidationError) as error:
        ln.Feature(
            name="feat_nested",
            dtype=list[ln.Record],
            cat_filters={"source__uid": "abc"},
        )
    assert (
        "lamindb.errors.ValidationError: cat_filters are incompatible with nested dtypes:"
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


def test_feature_from_dataframe():
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


def test_feature_from_dataframe_mute_restores_logger_verbosity():
    df = pd.DataFrame({"feature_mute_restore": [1, 2]})
    original_verbosity = logger._verbosity
    try:
        ln.Feature.from_dataframe(df, mute=True)
        assert logger._verbosity == original_verbosity
    finally:
        logger.set_verbosity(original_verbosity)


def test_feature_from_df_deprecation_warning():
    df = pd.DataFrame({"feature_from_df_deprecated": [1, 2]})
    with pytest.warns(DeprecationWarning, match="from_dataframe"):
        features = ln.Feature.from_df(df, mute=True)
    assert len(features) == 1


def test_feature_from_dict_mute_restores_logger_verbosity(dict_data):
    original_verbosity = logger._verbosity
    try:
        ln.Feature.from_dict(dict_data, mute=True)
        assert logger._verbosity == original_verbosity
    finally:
        logger.set_verbosity(original_verbosity)


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


def test_serialize_pandas_datetime_dtypes():
    datetime_series = pd.Series([pd.Timestamp("2024-01-01 12:00:00")])
    datetime_tz_series = pd.Series([pd.Timestamp("2024-01-01 12:00:00+00:00")])
    string_cat_series = pd.Series(["a", "b", "a"], dtype="category")

    assert serialize_pandas_dtype(datetime_series.dtype) == "datetime"
    assert serialize_pandas_dtype(datetime_tz_series.dtype) == "datetime64[ns, UTC]"
    assert serialize_pandas_dtype(string_cat_series.dtype) == "cat[ULabel]"


def test_dtype_as_object_covers_simple_fallbacks():
    assert dtype_as_object("datetime64[ns, UTC]") is datetime
    assert dtype_as_object("dict") is dict
    assert dtype_as_object("cat") is None
    assert dtype_as_object(None) is None  # type: ignore[arg-type]


def test_serialize_dtype_dict_and_invalid_type():
    assert serialize_dtype(dict) == "dict"
    with pytest.raises(
        ValueError,
        match="dtype has to be a registry, a ulabel subtype, a registry field",
    ):
        serialize_dtype(object())


def test_convert_to_pandas_dtype_unknown_roundtrip():
    assert convert_to_pandas_dtype("custom_dtype") == "custom_dtype"


def test_split_filter_parts_handles_escaped_commas():
    parts = _split_filter_parts(r"name='a\,b',status=active")
    assert parts == [r"name='a\,b'", "status=active"]


def test_format_cat_filter_value_edge_cases():
    assert _format_cat_filter_value(3.14) == "3.14"
    assert _format_cat_filter_value('a,"b') == "'a,\"b'"
    with pytest.raises(
        ValidationError,
        match="Cannot serialize categorical filter value containing comma and both quote types",
    ):
        _format_cat_filter_value("a,\"b'c")
