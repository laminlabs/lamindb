import re

import lamindb as ln
import pandas as pd
import pytest
from django.db import transaction
from lamindb.curators.core import DictCurator
from lamindb.errors import ValidationError


@pytest.fixture(autouse=True)
def isolate_records():
    with transaction.atomic():
        yield
        transaction.set_rollback(True)


def test_regex_roundtrip_and_reuse():
    feature = ln.Feature(name="regex_sample", dtype=str, regex="A+").save()
    assert ln.Feature.get(uid=feature.uid).regex == "A+"
    assert ln.Feature(name=feature.name, dtype=str).regex == "A+"
    assert ln.Feature(name=feature.name, dtype=str, regex=None).regex == "A+"
    replacement = ln.Feature(name=feature.name, dtype=str, regex="")
    assert replacement.regex == ""
    assert ln.Feature.get(uid=feature.uid).regex == "A+"
    replacement.save()
    assert ln.Feature.get(uid=feature.uid).regex == ""
    replacement._aux["unrelated"] = {"keep": True}
    replacement.regex = None
    replacement.save()
    restored = ln.Feature.get(uid=feature.uid)
    assert restored.regex is None
    assert restored._aux["unrelated"] == {"keep": True}


@pytest.mark.parametrize("pattern", [42, re.compile("A"), "["])
def test_regex_invalid_pattern_preserves_previous(pattern):
    feature = ln.Feature(name="regex_sample", dtype=str, regex="A")
    with pytest.raises((TypeError, ValidationError), match="regex"):
        feature.regex = pattern
    assert feature.regex == "A"


@pytest.mark.parametrize(
    "dtype", [int, float, bool, dict, list[str], "cat[ULabel]", "cat[ULabel|Record]"]
)
def test_regex_rejects_unsupported_dtype(dtype):
    with pytest.raises(ValidationError, match="regex"):
        ln.Feature(name="regex_sample", dtype=dtype, regex=".*")


def test_regex_definition_validation():
    with pytest.raises(ValidationError, match="regex"):
        ln.Feature(name="regex_type", is_type=True, regex="A")
    with pytest.raises(ValidationError, match="regex"):
        ln.Feature(
            name="regex_sample",
            dtype="cat[ULabel]",
            regex="A",
            cat_filters={"name": "A"},
        )
    for default in [1, "B"]:
        with pytest.raises(ValidationError, match="default_value"):
            ln.Feature(name="regex_sample", dtype=str, regex="A", default_value=default)
    feature = ln.Feature(
        name="regex_sample", dtype=str, regex="A", default_value="A"
    ).save()
    feature.regex = "B"
    with pytest.raises(ValidationError, match="default_value"):
        feature.save()
    assert ln.Feature.get(uid=feature.uid).regex == "A"
    feature.default_value = "B"
    feature.save()
    feature._dtype_str = "int"
    with pytest.raises(ValidationError, match="regex"):
        feature.save()
    assert ln.Feature.get(uid=feature.uid).dtype_as_str == "str"


@pytest.mark.parametrize("dtype", [str, "path", "url"])
@pytest.mark.parametrize(
    "value", ["sample(1)", "sample:1", "sample 1", "sample_1\n", ""]
)
def test_regex_fullmatch(dtype, value):
    feature = ln.Feature(
        name="regex_sample", dtype=dtype, regex=r"[A-Za-z0-9_-]+"
    ).save()
    schema = ln.Schema([feature]).save()
    ln.curators.DataFrameCurator(
        pd.DataFrame({feature.name: ["sample_1", None]}), schema
    ).validate()
    with pytest.raises(ValidationError, match="regex") as error:
        ln.curators.DataFrameCurator(
            pd.DataFrame({feature.name: ["sample_1", value]}, index=[10, 20]), schema
        ).validate()
    failures = error.value.__cause__.failure_cases
    assert failures["index"].tolist() == [20]
    assert failures["failure_case"].tolist() == [value]
    assert feature.name in str(error.value)


@pytest.mark.parametrize(
    "pattern,valid,invalid",
    [
        ("", "", "A"),
        (r"(?i)(sample)_\1", "Sample_SAMPLE", "Sample_other"),
        (r"(?m)^A$", "A", "A\n"),
    ],
)
def test_regex_python_semantics(pattern, valid, invalid):
    feature = ln.Feature(name="regex_sample", dtype=str, regex=pattern).save()
    schema = ln.Schema([feature]).save()
    DictCurator({feature.name: valid}, schema).validate()
    with pytest.raises(ValidationError, match="regex"):
        DictCurator({feature.name: invalid}, schema).validate()


def test_regex_index_and_annotation_updates():
    feature = ln.Feature(
        name="regex_sample", dtype=str, regex=r"[A-Za-z0-9_-]+", nullable=False
    ).save()
    schema = ln.Schema([], index=feature).save()
    ln.curators.DataFrameCurator(
        pd.DataFrame(index=pd.Index(["sample_1"], name=feature.name)), schema
    ).validate()
    with pytest.raises(ValidationError, match="regex"):
        ln.curators.DataFrameCurator(
            pd.DataFrame(index=pd.Index(["sample(1)"], name=feature.name)), schema
        ).validate()
    record_type = ln.Record(name="regex_sheet", is_type=True, schema=schema).save()
    record = ln.Record(name="sample_1", type=record_type).save()
    with pytest.raises(ValidationError, match="regex"):
        record.features.update_values({feature: "sample(1)"})
    record.refresh_from_db()
    assert record.name == "sample_1"
    record.name = "sample(1)"
    with pytest.raises(ValidationError, match="regex"):
        record.save()
    record.refresh_from_db()
    assert record.name == "sample_1"


@pytest.mark.parametrize(
    "field", ["name", "description", "reference", "reference_type"]
)
def test_regex_mapped_record_save(field):
    feature = ln.Feature(
        name="regex_sample", dtype=str, values_through=field, regex="A+"
    ).save()
    schema = ln.Schema([feature]).save()
    record_type = ln.Record(name="regex_sheet", is_type=True, schema=schema).save()
    record = ln.Record(type=record_type, **{field: "AA"}).save()
    setattr(record, field, "AB")
    with pytest.raises(ValidationError, match="regex"):
        record.save()
    record.refresh_from_db()
    assert getattr(record, field) == "AA"


def test_regex_annotation_without_schema():
    feature = ln.Feature(name="regex_sample", dtype=str, regex="A+").save()
    record = ln.Record(name="regex_row", features={feature: "AA"}).save()
    with pytest.raises(ValidationError, match="regex"):
        record.features.update_values({feature: "AB"})
    assert record.features.get_values()[feature.name] == "AA"


def test_regex_malformed_persisted_rule():
    feature = ln.Feature(name="regex_sample", dtype=str).save()
    schema = ln.Schema([feature]).save()
    ln.Feature.filter(uid=feature.uid).update(_aux={"regex": "["})
    assert ln.Feature.get(uid=feature.uid).regex == "["
    with pytest.raises(ValidationError, match="regex"):
        ln.curators.DataFrameCurator(pd.DataFrame({feature.name: ["A"]}), schema)


@pytest.mark.parametrize("lazy", [False, True])
def test_regex_polars(lazy):
    pl = pytest.importorskip("polars")
    feature = ln.Feature(
        name="regex_sample", dtype=str, regex=r"(?i)(sample)_\1"
    ).save()
    schema = ln.Schema([feature]).save()
    for value in ["Sample_SAMPLE", "Sample_other"]:
        df = pl.DataFrame({feature.name: [None, value]})
        curator = ln.curators.DataFrameCurator(df.lazy() if lazy else df, schema)
        if value == "Sample_SAMPLE":
            curator.validate()
        else:
            with pytest.raises(ValidationError, match="regex"):
                curator.validate()


@pytest.mark.parametrize("coerce", [None, False, True])
def test_regex_nulls_defaults_and_no_stringification(coerce):
    feature = ln.Feature(
        name="regex_sample",
        dtype=str,
        regex="A+",
        coerce=coerce,
        nullable=False,
        default_value="A",
    ).save()
    schema = ln.Schema([feature]).save()
    with pytest.raises(ValidationError):
        DictCurator({feature.name: None}, schema).validate()
    with pytest.raises(ValidationError):
        DictCurator({feature.name: 1}, schema).validate()
    curator = ln.curators.DataFrameCurator(pd.DataFrame({feature.name: [None]}), schema)
    curator.standardize()
    curator.validate()
    assert curator.dataset[feature.name].tolist() == ["A"]
    with pytest.raises(ValidationError, match="default_value"):
        ln.Feature(name=feature.name, dtype=str, regex="B")


@pytest.mark.parametrize(
    "field", ["created_at", "updated_at", "created_by", "run", "type", "unknown"]
)
def test_regex_incompatible_field_mapping(field):
    with pytest.raises((ValueError, ValidationError)):
        ln.Feature(name="regex_sample", dtype=str, regex="A", values_through=field)


def test_regex_bulk_record_values():
    feature = ln.Feature(name="regex_sample", dtype=str, regex="A+").save()
    schema = ln.Schema([feature]).save()
    record_type = ln.Record(name="regex_sheet", is_type=True, schema=schema).save()
    ln.Record.from_dataframe(
        pd.DataFrame({feature.name: ["A", "AA"]}), type=record_type
    ).save()
    with pytest.raises(ValidationError, match="regex"):
        ln.Record.from_dataframe(
            pd.DataFrame({feature.name: ["A", "AB"]}), type=record_type
        ).save()


def test_regex_curator_snapshot_and_removal():
    feature = ln.Feature(name="regex_sample", dtype=str, regex="A").save()
    schema = ln.Schema([feature]).save()
    data = pd.DataFrame({feature.name: ["A"]})
    old_curator = ln.curators.DataFrameCurator(data, schema)
    feature.regex = "B"
    feature.save()
    old_curator.validate()
    with pytest.raises(ValidationError, match="regex"):
        ln.curators.DataFrameCurator(data, schema).validate()
    feature.regex = None
    feature.save()
    ln.curators.DataFrameCurator(data, schema).validate()


def test_regex_artifact_and_run_annotations():
    feature = ln.Feature(name="regex_sample", dtype=str, regex="A+").save()
    transform = ln.Transform(key="regex-test-transform").save()
    run = ln.Run(transform).save()
    artifact = ln.Artifact.from_dataframe(
        pd.DataFrame({"value": [1]}), key="regex-test.parquet", run=run
    ).save()
    try:
        for host in [run, artifact]:
            host.features.set_values({feature: "AA"})
            with pytest.raises(ValidationError, match="regex"):
                host.features.update_values({feature: "AB"})
            assert host.features.get_values()[feature.name] == "AA"
    finally:
        artifact.delete(permanent=True, storage=True)


def test_regex_index_nullability_and_bulk_import():
    feature = ln.Feature(
        name="regex_sample", dtype=str, regex="A+", nullable=True
    ).save()
    schema = ln.Schema([], index=feature).save()
    frame = pd.DataFrame(index=pd.Index(["AA", None], name=feature.name))
    ln.curators.DataFrameCurator(frame, schema).validate()
    feature.nullable = False
    feature.save()
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(frame, schema).validate()
    record_type = ln.Record(name="regex_sheet", is_type=True, schema=schema).save()
    with pytest.raises(ValidationError, match="regex"):
        ln.Record.from_dataframe(
            pd.DataFrame(index=pd.Index(["AA", "AB"], name=feature.name)),
            type=record_type,
        ).save()


def test_regex_multiple_features_and_partial_record_save():
    first = ln.Feature(
        name="regex_first", dtype=str, regex="A+", values_through="name", nullable=False
    ).save()
    second = ln.Feature(
        name="regex_second", dtype=str, regex="B+", values_through="description"
    ).save()
    schema = ln.Schema([first, second]).save()
    DictCurator({first.name: "AA", second.name: "BB"}, schema).validate()
    with pytest.raises(ValidationError, match="regex_first"):
        DictCurator({first.name: "BB", second.name: "BB"}, schema).validate()
    record_type = ln.Record(name="regex_sheet", is_type=True, schema=schema).save()
    record = ln.Record(type=record_type, name="AA", description="BB").save()
    record.name = "invalid"
    record.description = "BBB"
    record.save(update_fields=["description"])
    record.refresh_from_db()
    assert record.name == "AA"
    assert record.description == "BBB"
    record.name = None
    with pytest.raises(ValidationError):
        record.save()


@pytest.mark.parametrize(
    "pattern,valid,invalid",
    [("", "", "A"), ("(?m)^A$", "A", "A\n"), (r"(?=A)A+", "AA", "AAB")],
)
def test_regex_polars_fullmatch(pattern, valid, invalid):
    pl = pytest.importorskip("polars")
    feature = ln.Feature(name="regex_sample", dtype=str, regex=pattern).save()
    schema = ln.Schema([feature]).save()
    ln.curators.DataFrameCurator(
        pl.DataFrame({feature.name: [valid, None]}), schema
    ).validate()
    with pytest.raises(ValidationError, match="regex"):
        ln.curators.DataFrameCurator(
            pl.DataFrame({feature.name: [invalid]}), schema
        ).validate()
