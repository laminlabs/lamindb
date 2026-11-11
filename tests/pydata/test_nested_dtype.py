from datetime import date, datetime

import lamindb as ln
import pandas as pd
import pytest
from django.db import transaction
from lamindb.errors import InvalidArgument, ValidationError
from lamindb.models.feature import (
    dtype_as_object,
    parse_dtype,
    serialize_dtype,
    transfer_feature_dtypes,
)


@pytest.fixture(autouse=True)
def isolate_records():
    with transaction.atomic():
        yield
        transaction.set_rollback(True)


@pytest.fixture(params=["pandas", "polars", "lazy"])
def make_frame(request):
    def frame(values):
        if request.param == "pandas":
            return pd.DataFrame({"nested": values})
        pl = pytest.importorskip("polars")
        result = pl.DataFrame({"nested": values}, strict=False)
        return result.lazy() if request.param == "lazy" else result

    return frame


def field_schema(*fields, **kwargs):
    return ln.Schema(features=[field.save() for field in fields], **kwargs).save()


def outer_schema(inner, *, dtype=None, nullable=True, **kwargs):
    return field_schema(
        ln.Feature(
            name="nested",
            dtype=inner if dtype is None else dtype,
            nullable=nullable,
        ),
        **kwargs,
    )


def test_dtype_serialization_and_reload(monkeypatch):
    inner = field_schema(ln.Feature(name="count", dtype=int))
    feature = ln.Feature(name="nested", dtype=inner).save()
    feature.refresh_from_db()
    assert feature.dtype_as_str == f"dict[{inner.uid}]"
    assert feature.dtype_as_object == inner
    assert serialize_dtype(list[list[inner]]) == f"list[list[dict[{inner.uid}]]]"
    assert dtype_as_object(serialize_dtype(list[inner])).__args__ == (inner,)
    assert parse_dtype(feature.dtype_as_str, check_exists=True) == []
    assert serialize_dtype(dict) == "dict"
    assert dtype_as_object("dict") is dict
    assert dtype_as_object("list[dict]") == list[dict]

    transfers = []
    monkeypatch.setattr(
        "lamindb.models.feature._transfer_dtype_schema",
        lambda uid, source, target, logs: transfers.append((uid, source, target)),
    )
    list_feature = ln.Feature(name="nested_list", dtype=list[inner]).save()
    transfer_feature_dtypes(list_feature, "target", {})
    assert transfers == [(inner.uid, list_feature._state.db, "target")]


def test_reject_unsupported_nested_schemas():
    unsaved = ln.Schema(features=[ln.Feature(name="count", dtype=int).save()])
    with pytest.raises(InvalidArgument, match="Save the nested Schema"):
        ln.Feature(name="nested", dtype=unsaved)
    flexible = ln.Schema(itype=ln.Feature).save()
    with pytest.raises(InvalidArgument, match="explicit Feature schema"):
        ln.Feature(name="nested", dtype=flexible)
    index = ln.Feature(name="index", dtype=int).save()
    indexed = ln.Schema(
        features=[ln.Feature(name="value", dtype=int).save()], index=index
    ).save()
    with pytest.raises(InvalidArgument, match="without slots or an index"):
        ln.Feature(name="nested", dtype=indexed)
    with pytest.raises(ValueError):
        parse_dtype("dict[not-a-schema]")
    with pytest.raises(ln.Schema.DoesNotExist):
        parse_dtype("list[dict[aAaAaAaAaAaAaAaA]]", check_exists=True)


def test_recursive_fields_and_optional_members(make_frame):
    details = field_schema(
        ln.Feature(name="count", dtype=int, nullable=False),
        ln.Feature(name="note", dtype=str),
        maximal_set=True,
    )
    note = details.members.get(name="note")
    details.optionals.add(note)
    details.save()
    inner = field_schema(
        ln.Feature(name="details", dtype=details, nullable=False),
        ln.Feature(name="tags", dtype=list[str]),
    )
    schema = outer_schema(inner)
    frame = make_frame([{"details": {"count": 2}, "tags": ["a"]}, None])
    curator = ln.curators.DataFrameCurator(frame, schema)
    curator.validate()
    assert type(curator.dataset) is type(frame)


@pytest.mark.parametrize(
    "values",
    [
        [{"count": "wrong"}],
        [{"other": 1}],
        [{}],
        [{"count": None}],
        [{"count": 1, "other": 2}],
        ["not a dictionary"],
    ],
)
def test_invalid_fields(make_frame, values):
    inner = field_schema(
        ln.Feature(name="count", dtype=int, nullable=False), maximal_set=True
    )
    curator = ln.curators.DataFrameCurator(make_frame(values), outer_schema(inner))
    with pytest.raises(ValidationError, match="nested"):
        curator.validate()
    assert not curator._is_validated


def test_required_nullable_field_uses_shared_layout(make_frame):
    inner = field_schema(ln.Feature(name="count", dtype=int))
    curator = ln.curators.DataFrameCurator(
        make_frame([{"count": 1}, {}]), outer_schema(inner)
    )
    curator.validate()
    ln.curators.DataFrameCurator(
        make_frame([{"count": None}]), outer_schema(inner)
    ).validate()
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(make_frame([{}]), outer_schema(inner)).validate()


@pytest.mark.parametrize("nullable", [False, True])
def test_parent_nullability(make_frame, nullable):
    inner = field_schema(ln.Feature(name="count", dtype=int, nullable=False))
    curator = ln.curators.DataFrameCurator(
        make_frame([{"count": 1}, None]), outer_schema(inner, nullable=nullable)
    )
    if nullable:
        curator.validate()
    else:
        with pytest.raises(ValidationError):
            curator.validate()


@pytest.mark.parametrize("depth", [1, 2])
def test_list_composition(make_frame, depth):
    inner = field_schema(ln.Feature(name="count", dtype=int, nullable=False))
    dtype = list[inner] if depth == 1 else list[list[inner]]
    value = [{"count": 1}, None] if depth == 1 else [[{"count": 1}], []]
    schema = outer_schema(inner, dtype=dtype)
    ln.curators.DataFrameCurator(make_frame([value, [], None]), schema).validate()
    wrong = [{"count": "wrong"}] if depth == 1 else [[{"count": "wrong"}]]
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(make_frame([wrong]), schema).validate()


def test_nested_lists_of_scalars(make_frame):
    inner = field_schema(ln.Feature(name="counts", dtype=list[int]))
    ln.curators.DataFrameCurator(
        make_frame([{"counts": [1, 2]}]), outer_schema(inner)
    ).validate()
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            make_frame([{"counts": ["wrong"]}]), outer_schema(inner)
        ).validate()


def test_nested_categorical_validation(make_frame):
    ln.ULabel(name="nested known label").save()
    inner = field_schema(ln.Feature(name="label", dtype=ln.ULabel))
    schema = outer_schema(inner)
    ln.curators.DataFrameCurator(
        make_frame([{"label": "nested known label"}]), schema
    ).validate()
    with pytest.raises(ValidationError, match="nested"):
        ln.curators.DataFrameCurator(
            make_frame([{"label": "nested unknown label"}]), schema
        ).validate()


def test_nested_minimal_set_false(make_frame):
    inner = field_schema(ln.Feature(name="count", dtype=int), minimal_set=False)
    ln.curators.DataFrameCurator(
        make_frame([{"other": "allowed"}]), outer_schema(inner)
    ).validate()


@pytest.mark.parametrize("coerce_on", ["member", "feature", "inner", "outer"])
def test_nested_coercion_rejected(make_frame, coerce_on):
    inner = field_schema(
        ln.Feature(name="count", dtype=int, coerce=coerce_on == "member"),
        coerce=coerce_on == "inner",
    )
    schema = field_schema(
        ln.Feature(name="nested", dtype=inner, coerce=coerce_on == "feature"),
        coerce=coerce_on == "outer",
    )
    with pytest.raises(InvalidArgument, match="coerce|coercion"):
        ln.curators.DataFrameCurator(make_frame([{"count": "1"}]), schema).validate()


def test_nested_scalar_types(make_frame):
    inner = field_schema(
        ln.Feature(name="count", dtype=int),
        ln.Feature(name="fraction", dtype=float),
        ln.Feature(name="flag", dtype=bool),
        ln.Feature(name="day", dtype=date),
        ln.Feature(name="timestamp", dtype=datetime),
        ln.Feature(name="measurement", dtype="num"),
    )
    value = {
        "count": 1,
        "fraction": 0.5,
        "flag": True,
        "day": date(2026, 1, 1),
        "timestamp": datetime(2026, 1, 1),
        "measurement": 2,
    }
    ln.curators.DataFrameCurator(make_frame([value]), outer_schema(inner)).validate()


@pytest.mark.parametrize("nullable", [False, True])
def test_nested_nan_nullability(make_frame, nullable):
    inner = field_schema(ln.Feature(name="fraction", dtype=float, nullable=nullable))
    curator = ln.curators.DataFrameCurator(
        make_frame([{"fraction": float("nan")}]), outer_schema(inner)
    )
    if nullable:
        curator.validate()
    else:
        with pytest.raises(ValidationError):
            curator.validate()


def test_optional_nonnullable_field(make_frame):
    inner = field_schema(
        ln.Feature(name="count", dtype=int),
        ln.Feature(name="note", dtype=str, nullable=False),
    )
    inner.optionals.add(inner.members.get(name="note"))
    inner.save()
    schema = outer_schema(inner)
    ln.curators.DataFrameCurator(make_frame([{"count": 1}]), schema).validate()
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            make_frame([{"count": 1, "note": None}]), schema
        ).validate()


def test_nested_field_order(make_frame):
    inner = field_schema(
        ln.Feature(name="count", dtype=int),
        ln.Feature(name="note", dtype=str),
        ordered_set=True,
    )
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            make_frame([{"note": "wrong order", "count": 1}]), outer_schema(inner)
        ).validate()


def test_optional_nonnullable_field_uses_shared_layout(make_frame):
    inner = field_schema(
        ln.Feature(name="count", dtype=int),
        ln.Feature(name="note", dtype=str, nullable=False),
    )
    inner.optionals.add(inner.members.get(name="note"))
    inner.save()
    frame = make_frame([{"count": 1}, {"count": 2, "note": "present"}])
    curator = ln.curators.DataFrameCurator(frame, outer_schema(inner))
    with pytest.raises(ValidationError):
        curator.validate()


@pytest.mark.parametrize(
    "values,options,valid",
    [
        ([{}, {"count": 1}], {"minimal_set": True}, True),
        ([{"note": "a"}, {"note": "b"}], {"minimal_set": True}, False),
        ([{"note": "a"}, {"note": "b"}], {"minimal_set": False}, True),
        ([{}, {"count": 1, "extra": 2}], {"maximal_set": True}, False),
        ([{}, {"count": 1, "extra": 2}], {"maximal_set": False}, True),
        (
            [{"count": 1, "note": "a"}, {"note": "b", "count": 2}],
            {"ordered_set": True},
            True,
        ),
        (
            [{"note": "a", "count": 1}, {"count": 2, "note": "b"}],
            {"ordered_set": True},
            False,
        ),
        (
            [{"note": "a", "count": 1}, {"count": 2, "note": "b"}],
            {"ordered_set": False},
            True,
        ),
        (
            [{"note": "a"}, {"count": 2, "note": "b"}],
            {"ordered_set": True},
            False,
        ),
        ([{"count": 1}, {"note": "a"}], {"ordered_set": True}, True),
    ],
)
def test_shared_layout_set_constraints(make_frame, values, options, valid):
    inner = field_schema(
        ln.Feature(name="count", dtype=int),
        ln.Feature(name="note", dtype=str),
        **options,
    )
    inner.optionals.add(inner.members.get(name="note"))
    inner.save()
    frame = make_frame(values)
    if "ordered_set" in options and not isinstance(frame, pd.DataFrame):
        import polars as pl

        # Compare the same layout, not Polars' constructor-dependent inference order.
        names = dict.fromkeys(name for value in values for name in value)
        frame = frame.with_columns(
            pl.struct([pl.col("nested").struct.field(name) for name in names]).alias(
                "nested"
            )
        )
    curator = ln.curators.DataFrameCurator(frame, outer_schema(inner))
    if valid:
        curator.validate()
    else:
        with pytest.raises(ValidationError):
            curator.validate()


@pytest.mark.parametrize("as_list", [False, True])
@pytest.mark.parametrize("nullable", [False, True])
def test_recursive_shared_layout(make_frame, as_list, nullable):
    details = field_schema(
        ln.Feature(name="count", dtype=int),
        ln.Feature(name="note", dtype=str, nullable=nullable),
    )
    details.optionals.add(details.members.get(name="note"))
    details.save()
    inner = field_schema(
        ln.Feature(name="details", dtype=list[details] if as_list else details)
    )
    values = [
        {"details": [{"count": 1}]} if as_list else {"details": {"count": 1}},
        None,
        (
            {"details": [{"note": "a", "count": 2}]}
            if as_list
            else {"details": {"note": "a", "count": 2}}
        ),
    ]
    curator = ln.curators.DataFrameCurator(make_frame(values), outer_schema(inner))
    if nullable:
        curator.validate()
    else:
        with pytest.raises(ValidationError):
            curator.validate()


@pytest.mark.parametrize("depth", [1, 2])
@pytest.mark.parametrize("nullable", [False, True])
def test_shared_layout_across_list_elements(make_frame, depth, nullable):
    inner = field_schema(
        ln.Feature(name="count", dtype=int),
        ln.Feature(name="note", dtype=str, nullable=nullable),
    )
    inner.optionals.add(inner.members.get(name="note"))
    inner.save()
    first = [{}, {"count": 1}]
    second = [{"count": 2, "note": "present"}]
    values = (
        [first, [], None, second] if depth == 1 else [[first, []], [], None, [second]]
    )
    dtype = list[inner] if depth == 1 else list[list[inner]]
    curator = ln.curators.DataFrameCurator(
        make_frame(values), outer_schema(inner, dtype=dtype)
    )
    if nullable:
        curator.validate()
    else:
        with pytest.raises(ValidationError):
            curator.validate()


@pytest.mark.parametrize(
    "dtype,value",
    [
        (int, 1),
        (float, 1.5),
        (bool, True),
        (str, "value"),
        (date, date(2026, 1, 1)),
        (datetime, datetime(2026, 1, 1)),
        (list[int], [1]),
        (dict, {"unconstrained": 1}),
    ],
)
def test_missing_nullable_scalar_field(make_frame, dtype, value):
    inner = field_schema(ln.Feature(name="value", dtype=dtype))
    ln.curators.DataFrameCurator(
        make_frame([{}, {"value": value}, {"value": None}]), outer_schema(inner)
    ).validate()


def test_shared_layout_does_not_mutate_source():
    inner = field_schema(
        ln.Feature(name="count", dtype=int),
        ln.Feature(name="note", dtype=str),
        ordered_set=True,
    )
    values = [{"count": 1, "note": "a"}, {"note": "b", "count": 2}, {}]
    frame = pd.DataFrame({"nested": values})
    before = frame.copy(deep=True)
    curator = ln.curators.DataFrameCurator(frame, outer_schema(inner))
    curator.validate()
    pd.testing.assert_frame_equal(curator.dataset, before)
    assert list(values[1]) == ["note", "count"]
    assert values[2] == {}


def test_nested_bare_dict_string_compatibility():
    inner = field_schema(ln.Feature(name="value", dtype=dict))
    ln.curators.DataFrameCurator(
        pd.DataFrame({"nested": [{}, {"value": "{'unconstrained': 1}"}]}),
        outer_schema(inner),
    ).validate()


@pytest.mark.parametrize("minimal_set", [False, True])
@pytest.mark.parametrize("optional", [False, True])
@pytest.mark.parametrize("nullable", [False, True])
def test_field_presence_and_nullability_are_independent(
    make_frame, minimal_set, optional, nullable
):
    value = ln.Feature(name="value", dtype=int, nullable=nullable).save()
    inner = field_schema(
        ln.Feature(name="count", dtype=int),
        value,
        minimal_set=minimal_set,
    )
    if optional:
        inner.optionals.add(value)
        inner.save()
    schema = outer_schema(inner)
    mixed = ln.curators.DataFrameCurator(
        make_frame([{"count": 1}, {"count": 2, "value": 3}]), schema
    )
    absent = ln.curators.DataFrameCurator(make_frame([{"count": 1}]), schema)
    if nullable:
        mixed.validate()
    else:
        with pytest.raises(ValidationError):
            mixed.validate()
    if minimal_set and not optional:
        with pytest.raises(ValidationError):
            absent.validate()
    else:
        absent.validate()


def test_optional_nested_column(make_frame):
    inner = field_schema(ln.Feature(name="count", dtype=int))
    schema = outer_schema(inner)
    schema.optionals.add(schema.members.get(name="nested"))
    schema.save()
    frame = make_frame([{"count": 1}])
    frame = (
        frame.drop(columns="nested")
        if isinstance(frame, pd.DataFrame)
        else frame.drop("nested")
    )
    ln.curators.DataFrameCurator(frame, schema).validate()


def test_all_null_parents(make_frame):
    inner = field_schema(ln.Feature(name="count", dtype=int, nullable=False))
    ln.curators.DataFrameCurator(
        make_frame([None, None]), outer_schema(inner)
    ).validate()


def test_empty_column(make_frame):
    inner = field_schema(ln.Feature(name="count", dtype=int, nullable=False))
    ln.curators.DataFrameCurator(
        make_frame([]), outer_schema(inner, nullable=False)
    ).validate()


def test_empty_nested_schema(make_frame):
    inner = ln.Schema(itype=ln.Feature, flexible=False, maximal_set=True).save()
    schema = outer_schema(inner)
    ln.curators.DataFrameCurator(make_frame([{}]), schema).validate()
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(make_frame([{"unexpected": 1}]), schema).validate()


def test_scoped_feature_schema(make_frame):
    feature_type = ln.Feature(name="Nested measurements", is_type=True).save()
    count = ln.Feature(name="count", dtype=int, type=feature_type).save()
    inner = ln.Schema(features=[count], itype=feature_type).save()
    ln.curators.DataFrameCurator(
        make_frame([{"count": 1}]), outer_schema(inner)
    ).validate()


def test_cycle_rejected():
    member = ln.Feature(name="count", dtype=int).save()
    inner = field_schema(member)
    ln.Feature.filter(uid=member.uid).update(_dtype_str=f"dict[{inner.uid}]")
    with pytest.raises(InvalidArgument, match="Cyclic nested dtype"):
        ln.curators.DataFrameCurator(
            pd.DataFrame({"nested": [{"count": {}}]}), outer_schema(inner)
        ).validate()


@pytest.mark.parametrize("as_list", [False, True])
def test_record_values_roundtrip(as_list):
    inner = field_schema(ln.Feature(name="count", dtype=int, nullable=False))
    dtype = list[inner] if as_list else inner
    feature = ln.Feature(name="nested", dtype=dtype).save()
    value = [{"count": 2}] if as_list else {"count": 2}
    record = ln.Record(name="nested record", features={feature: value}).save()
    record.refresh_from_db()
    assert record.features.get_values()["nested"] == value
    wrong = [{"count": "wrong"}] if as_list else {"count": "wrong"}
    with pytest.raises(ValidationError):
        record.features.set_values({feature: wrong})
    assert record.features.get_values()["nested"] == value


@pytest.mark.parametrize("as_list", [False, True])
def test_curated_artifact_roundtrip(make_frame, as_list):
    inner = field_schema(ln.Feature(name="count", dtype=int, nullable=False))
    schema = outer_schema(inner, dtype=list[inner] if as_list else inner)
    value = [{"count": 1}] if as_list else {"count": 1}
    curator = ln.curators.DataFrameCurator(make_frame([value]), schema)
    artifact = curator.save_artifact(key="nested/curated.parquet")
    artifact.refresh_from_db()
    assert artifact.schema == schema
    loaded = artifact.load()
    ln.curators.DataFrameCurator(loaded, schema).validate()


def test_revalidation_resets_state():
    inner = field_schema(ln.Feature(name="count", dtype=int))
    frame = pd.DataFrame({"nested": [{"count": 1}]})
    curator = ln.curators.DataFrameCurator(frame, outer_schema(inner))
    curator.validate()
    frame.loc[0, "nested"]["count"] = "wrong"
    with pytest.raises(ValidationError):
        curator.validate()
    assert not curator._is_validated


@pytest.mark.parametrize("invalid", [None, "parent-null", "child-null", "wrong-type"])
def test_polars_lazy_projection_only(tmp_path, monkeypatch, invalid):
    pl = pytest.importorskip("polars")
    inner = field_schema(ln.Feature(name="count", dtype=int, nullable=False))
    schema = outer_schema(inner, nullable=False)
    path = tmp_path / "nested.parquet"
    values = [{"count": i} for i in range(1000)]
    if invalid == "parent-null":
        values[0] = None
    elif invalid == "child-null":
        values[0] = {"count": None}
    elif invalid == "wrong-type":
        values = [{"count": "wrong"} for _ in range(1000)]
    pl.DataFrame({"nested": values}).write_parquet(path)
    original_collect = pl.LazyFrame.collect
    collected_shapes = []

    def collect(frame, *args, **kwargs):
        result = original_collect(frame, *args, **kwargs)
        collected_shapes.append(result.shape)
        assert result.height <= 1, "Validation materialized dataset rows"
        return result

    monkeypatch.setattr(pl.LazyFrame, "collect", collect)
    curator = ln.curators.DataFrameCurator(pl.scan_parquet(path), schema)
    if invalid is None:
        curator.validate()
    else:
        with pytest.raises(ValidationError):
            curator.validate()
    assert isinstance(curator.dataset, pl.LazyFrame)
    assert collected_shapes


def test_bare_dict_compatibility(make_frame):
    schema = outer_schema(None, dtype=dict)
    ln.curators.DataFrameCurator(
        make_frame([{"anything": {"not": ["constrained"]}}]), schema
    ).validate()
