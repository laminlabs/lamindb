from datetime import date, datetime

import bionty as bt
import lamindb as ln
import pytest
from django.db import transaction
from lamindb.errors import InvalidArgument, ValidationError
from upath import UPath

pl = pytest.importorskip("polars")
pa = pytest.importorskip("pandera.polars")


@pytest.fixture(autouse=True)
def isolate_records():
    with transaction.atomic():
        yield
        transaction.set_rollback(True)


@pytest.fixture(params=[False, True], ids=["dataframe", "lazyframe"])
def as_frame(request):
    return lambda df: df.lazy() if request.param else df


def make_schema(features, **kwargs):
    return ln.Schema(
        features=[
            ln.Feature(name=name, dtype=dtype, **options).save()
            for name, dtype, options in features
        ],
        **kwargs,
    ).save()


def eager(frame):
    return frame.collect() if isinstance(frame, pl.LazyFrame) else frame


def test_polars_scalar_and_nested_dtypes(as_frame):
    schema = make_schema(
        [
            ("pl_int", int, {}),
            ("pl_float", float, {}),
            ("pl_str", str, {}),
            ("pl_bool", bool, {}),
            ("pl_num", "num", {}),
            ("pl_list", list[str], {}),
            ("pl_list_int", list[int], {}),
            ("pl_list_float", list[float], {}),
            ("pl_dict", dict, {}),
            ("pl_date", date, {}),
            ("pl_datetime", datetime, {}),
        ],
        maximal_set=True,
        ordered_set=True,
    )
    df = pl.DataFrame(
        {
            "pl_int": pl.Series([1, 2], dtype=pl.Int16),
            "pl_float": pl.Series([1.0, 2.0], dtype=pl.Float32),
            "pl_str": ["a", "b"],
            "pl_bool": [True, False],
            "pl_num": [1, 2],
            "pl_list": [["a"], ["b"]],
            "pl_list_int": pl.Series([[1], [2]], dtype=pl.List(pl.Int16)),
            "pl_list_float": pl.Series([[1.0], [2.0]], dtype=pl.List(pl.Float32)),
            "pl_dict": [{"a": 1}, {"a": 2}],
            "pl_date": [date(2024, 1, 1)] * 2,
            "pl_datetime": [datetime(2024, 1, 1)] * 2,
        }
    )
    curator = ln.curators.DataFrameCurator(as_frame(df), schema)
    assert type(curator.dataset) is type(as_frame(df))
    assert isinstance(curator._atomic_curator._pandera_schema, pa.DataFrameSchema)
    curator.validate()
    assert type(curator.dataset) is type(as_frame(df))
    assert eager(curator.dataset).equals(df)
    assert curator.dataset.collect_schema()["pl_int"] == pl.Int16
    assert curator.dataset.collect_schema()["pl_float"] == pl.Float32
    assert curator.dataset.collect_schema()["pl_list_int"] == pl.List(pl.Int16)
    assert curator.dataset.collect_schema()["pl_list_float"] == pl.List(pl.Float32)


@pytest.mark.parametrize(
    "dtype,values",
    [
        (int, ["a"]),
        (float, [1]),
        (str, [1]),
        (bool, [1]),
        (list[str], [[1]]),
        (list[int], [["a"]]),
        (list[float], [["a"]]),
        (dict, ["a"]),
        ("list[object]", [["a"]]),
    ],
)
@pytest.mark.parametrize("regex", [None, "a+"])
def test_polars_invalid_dtype(as_frame, dtype, values, regex):
    features = [("pl_invalid", dtype, {"regex": regex} if dtype is str else {})]
    data = {"pl_invalid": values}
    if regex is not None:
        features.append(("pl_regex", str, {"regex": regex}))
        data["pl_regex"] = ["b"]
    schema = make_schema(features)
    curator = ln.curators.DataFrameCurator(as_frame(pl.DataFrame(data)), schema)
    with pytest.raises(ValidationError) as error:
        curator.validate()
    assert "must fully match regex" not in str(error.value)
    assert not curator._is_validated


@pytest.mark.parametrize("coerce_on", ["feature", "schema"])
@pytest.mark.parametrize(
    "dtype,values,expected_dtype,expected_values",
    [
        (int, ["1", "2"], pl.Int64, [1, 2]),
        (float, ["1", "2"], pl.Float64, [1.0, 2.0]),
        (float, pl.Series([1.0, 2.0], dtype=pl.Float32), pl.Float32, [1.0, 2.0]),
    ],
    ids=["int-from-str", "float-from-str", "float32-preserved"],
)
def test_polars_coercion(
    as_frame, coerce_on, dtype, values, expected_dtype, expected_values
):
    schema = make_schema(
        [("pl_coerce", dtype, {"coerce": coerce_on == "feature"})],
        coerce=coerce_on == "schema",
    )
    df = pl.DataFrame({"pl_coerce": values})
    curator = ln.curators.DataFrameCurator(as_frame(df), schema)
    curator.validate()
    assert type(curator.dataset) is type(as_frame(df))
    assert curator.dataset.collect_schema()["pl_coerce"] == expected_dtype
    assert eager(curator.dataset)["pl_coerce"].to_list() == expected_values


def test_polars_coercion_is_lossless(as_frame):
    schema = make_schema([("pl_lossless", int, {"coerce": True})])
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_lossless": [1.5]})), schema
    )
    with pytest.raises(ValidationError):
        curator.validate()


@pytest.mark.parametrize("nullable", [True, False])
@pytest.mark.parametrize("value", ["a", "b"])
def test_polars_nullability(as_frame, nullable, value):
    schema = make_schema([("pl_nullable", str, {"nullable": nullable, "regex": "a+"})])
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_nullable": [value, None]})), schema
    )
    if nullable and value == "a":
        curator.validate()
    else:
        with pytest.raises(
            ValidationError, match="regex" if nullable else "non-nullable"
        ):
            curator.validate()


def test_polars_column_constraints(as_frame):
    schema = make_schema(
        [("pl_first", int, {}), ("pl_second", str, {"regex": "a+"})],
        maximal_set=True,
        ordered_set=True,
    )
    for df in [
        pl.DataFrame({"pl_first": [1]}),
        pl.DataFrame({"pl_second": ["a"], "pl_first": [1]}),
        pl.DataFrame({"pl_first": [1], "pl_second": ["a"], "extra": [1]}),
        pl.DataFrame({"pl_first": [1], "pl_second": ["ab"]}),
    ]:
        with pytest.raises(ValidationError):
            ln.curators.DataFrameCurator(as_frame(df), schema).validate()
    ln.curators.DataFrameCurator(
        as_frame(
            pl.DataFrame(
                {"pl_first": [1], "pl_second": ["a"], "__lamindb_record_uid__": ["x"]}
            )
        ),
        schema,
    ).validate()


def test_polars_optional_column(as_frame):
    required = ln.Feature(name="pl_required", dtype=int).save()
    optional = ln.Feature(name="pl_optional", dtype=str, regex="a+").save()
    schema = ln.Schema(features=[required, optional]).save()
    schema.optionals.add(optional)
    ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_required": [1]})), schema
    ).validate()


@pytest.mark.parametrize("dtype", [pl.String, pl.Categorical])
def test_polars_registry_validation(as_frame, dtype):
    ln.ULabel(name="polars known label").save()
    schema = make_schema([("pl_label", ln.ULabel, {})])
    df = pl.DataFrame(
        {"pl_label": pl.Series(["polars known label", "polars new label"], dtype=dtype)}
    )
    curator = ln.curators.DataFrameCurator(as_frame(df), schema)
    with pytest.raises(ValidationError):
        curator.validate()
    assert curator.cat.non_validated == {"pl_label": ["polars new label"]}
    curator.cat.add_new_from("pl_label")
    curator.validate()
    assert curator.cat.non_validated == {}


def test_polars_list_registry_validation(as_frame):
    ln.ULabel(name="polars list label").save()
    schema = make_schema([("pl_labels", list[ln.ULabel], {})])
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_labels": [["polars list label"]]})), schema
    )
    curator.validate()


def test_polars_standardize(as_frame):
    schema = make_schema(
        [
            ("pl_default", str, {"default_value": "filled"}),
            ("pl_missing", int, {"nullable": True}),
        ]
    )
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_default": ["a", None]})), schema
    )
    curator.standardize()
    assert type(curator.dataset) is type(as_frame(pl.DataFrame()))
    assert eager(curator.dataset)["pl_default"].to_list() == ["a", "filled"]
    assert eager(curator.dataset)["pl_missing"].to_list() == [None, None]
    curator.validate()

    required = make_schema([("pl_required_missing", int, {"nullable": False})])
    with pytest.raises(ValidationError, match="pl_required_missing"):
        ln.curators.DataFrameCurator(
            as_frame(pl.DataFrame({"other": [1]})), required
        ).standardize()


@pytest.mark.parametrize(
    "dtype,expected",
    [
        (date, pl.Date),
        (datetime, pl.Datetime()),
        ("datetime64[ns, UTC]", pl.Datetime(time_unit="ns", time_zone="UTC")),
    ],
)
def test_polars_standardize_nullable_temporal_column(as_frame, dtype, expected):
    schema = make_schema([("pl_missing_temporal", dtype, {"nullable": True})])
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_other": [1]})), schema
    )
    curator.standardize()
    assert curator.dataset.collect_schema()["pl_missing_temporal"] == expected
    assert eager(curator.dataset)["pl_missing_temporal"].to_list() == [None]
    curator.validate()


def test_polars_index_not_supported(as_frame):
    index = ln.Feature(name="pl_index", dtype=str).save()
    schema = ln.Schema(itype=ln.Feature, index=index).save()
    with pytest.raises(InvalidArgument, match="do not have an index"):
        ln.curators.DataFrameCurator(as_frame(pl.DataFrame({"a": [1]})), schema)


def test_polars_save_artifact(as_frame):
    schema = make_schema([("pl_saved", int, {})])
    df = pl.DataFrame({"pl_saved": [1, 2]})
    curator = ln.curators.DataFrameCurator(as_frame(df), schema)
    artifact = curator.save_artifact(key="polars/curated.parquet")
    assert artifact.schema == schema
    assert artifact.n_observations == 2
    assert eager(artifact.load())["pl_saved"].to_list() == [1, 2]
    assert ln.Artifact.get(uid=artifact.uid).load()["pl_saved"].tolist() == [1, 2]
    with artifact.open(engine="polars") as lazy:
        ln.curators.DataFrameCurator(lazy, schema).validate()


def test_polars_artifact_constructor(as_frame):
    schema = make_schema([("pl_constructor", int, {})])
    df = pl.DataFrame({"pl_constructor": [1, 2]})
    artifact = ln.Artifact.from_dataframe(
        as_frame(df), key="polars/constructor.parquet", schema=schema
    )
    assert isinstance(
        artifact._curator._atomic_curator._pandera_schema, pa.DataFrameSchema
    )
    artifact.save()
    assert artifact.n_observations == 2
    assert eager(artifact.load())["pl_constructor"].to_list() == [1, 2]
    with pytest.raises(ValidationError):
        ln.Artifact.from_dataframe(
            as_frame(pl.DataFrame({"pl_constructor": ["invalid"]})),
            key="polars/invalid.parquet",
            schema=schema,
        )


def test_lazyframe_construction_does_not_execute(monkeypatch, tmp_path):
    schema = make_schema([("pl_scanned", int, {})])
    path = tmp_path / "lazy.parquet"
    pl.DataFrame({"pl_scanned": [1, 2]}).write_parquet(path)
    lazy = pl.scan_parquet(path)

    def unexpected_collect(*args, **kwargs):
        pytest.fail("Constructing the curator must not execute the LazyFrame")

    with monkeypatch.context() as context:
        context.setattr(pl.LazyFrame, "collect", unexpected_collect)
        curator = ln.curators.DataFrameCurator(lazy, schema)
        assert curator.dataset is lazy
    curator.validate()
    assert isinstance(curator.dataset, pl.LazyFrame)
    assert "Parquet SCAN" in curator.dataset.explain()


@pytest.mark.parametrize("dtype", [int, float, str, bool, list[str]])
@pytest.mark.parametrize("values", [[], [None, None]])
def test_polars_empty_nullable_columns(as_frame, dtype, values):
    schema = make_schema([("pl_empty", dtype, {"nullable": True})])
    ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_empty": values})), schema
    ).validate()


def test_polars_empty_lists(as_frame):
    schema = make_schema([("pl_empty_list", list[str], {})])
    ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_empty_list": [[], []]})), schema
    ).validate()


@pytest.mark.parametrize("dtype", [int, float])
def test_polars_invalid_coercion(as_frame, dtype):
    schema = make_schema([("pl_bad_cast", dtype, {"coerce": True})])
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            as_frame(pl.DataFrame({"pl_bad_cast": ["not a number"]})), schema
        ).validate()


def test_polars_standardize_column_synonym(as_frame):
    feature = ln.Feature(name="pl_canonical_col", dtype=int).save()
    feature.add_synonym("pl_alias_col")
    schema = ln.Schema(features=[feature]).save()
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_alias_col": [1]})), schema
    )
    curator.cat.standardize("columns")
    assert list(curator.dataset.collect_schema().names()) == ["pl_canonical_col"]


def test_polars_biorecord_source_filter():
    from lamindb.curators.core import DataFrameCatManager

    source = bt.Source(
        entity="bionty.CellType", organism="all", name="polars-cl", version="1"
    ).save()
    schema = make_schema([("pl_cell_source", bt.CellType, {})])
    manager = DataFrameCatManager(
        pl.DataFrame({"pl_cell_source": ["T cell"]}),
        categoricals=list(schema.members),
        sources={"pl_cell_source": source},
        schema=schema,
    )
    assert manager._cat_vectors["pl_cell_source"]._filter_kwargs["source"] == source


def test_polars_unique_values_without_registry_columns():
    schema = make_schema([("pl_no_cat", int, {})])
    curator = ln.curators.DataFrameCurator(pl.DataFrame({"pl_no_cat": [1]}), schema)
    assert curator.cat._get_unique_values() == {}


@pytest.mark.parametrize("list_values", [False, True])
def test_polars_registry_synonyms(as_frame, list_values):
    label = bt.CellType(name="polars canonical", synonyms="polars synonym").save()
    dtype = list[bt.CellType] if list_values else bt.CellType
    schema = make_schema([("pl_synonym", dtype, {})])
    values = [["polars synonym"]] if list_values else ["polars synonym"]
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_synonym": values})), schema
    )
    with pytest.raises(ValidationError):
        curator.validate()
    curator.cat.standardize("pl_synonym")
    assert type(curator.dataset) is type(as_frame(pl.DataFrame()))
    expected = [[label.name]] if list_values else [label.name]
    assert eager(curator.dataset)["pl_synonym"].to_list() == expected
    curator.validate()


@pytest.mark.parametrize("list_values", [False, True])
def test_polars_enum_synonyms(as_frame, list_values):
    label = bt.CellType(
        name="polars enum canonical", synonyms="polars enum synonym"
    ).save()
    dtype = list[bt.CellType] if list_values else bt.CellType
    schema = make_schema([("pl_enum_synonym", dtype, {})])
    values = ["polars enum synonym"]
    enum = pl.Enum(values)
    series_dtype = pl.List(enum) if list_values else enum
    series_values = [["polars enum synonym"]] if list_values else values
    curator = ln.curators.DataFrameCurator(
        as_frame(
            pl.DataFrame(
                {"pl_enum_synonym": pl.Series(series_values, dtype=series_dtype)}
            )
        ),
        schema,
    )
    with pytest.raises(ValidationError):
        curator.validate()
    curator.cat.standardize("pl_enum_synonym")
    expected = [[label.name]] if list_values else [label.name]
    assert eager(curator.dataset)["pl_enum_synonym"].to_list() == expected
    restored = curator.dataset.collect_schema()["pl_enum_synonym"]
    restored_enum = restored.inner if list_values else restored
    assert isinstance(restored_enum, pl.Enum)
    assert list(restored_enum.categories) == ["polars enum synonym", label.name]
    curator.validate()


def test_polars_standardize_missing_registry_column(as_frame):
    schema = make_schema(
        [
            ("pl_present", int, {}),
            (
                "pl_missing_label",
                ln.ULabel,
                {"default_value": "unregistered polars label"},
            ),
        ]
    )
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_present": [1]})), schema
    )
    curator.standardize()
    with pytest.raises(ValidationError, match="unregistered polars label"):
        curator.validate()


def test_polars_revalidate_after_filling_unregistered_default(as_frame):
    ln.ULabel(name="polars cached label").save()
    schema = make_schema(
        [
            (
                "pl_cached_label",
                ln.ULabel,
                {
                    "nullable": True,
                    "default_value": "unregistered cached default",
                },
            )
        ]
    )
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_cached_label": ["polars cached label", None]})),
        schema,
    )
    curator.validate()
    assert curator.cat._unique_values is not None
    curator.standardize()
    assert curator.cat._unique_values is None
    with pytest.raises(ValidationError, match="unregistered cached default"):
        curator.validate()
    assert not ln.ULabel.filter(name="unregistered cached default").exists()


def test_polars_revalidate_after_adding_categorical_column(as_frame):
    known = ln.ULabel(name="polars present label").save()
    present = ln.Feature(name="pl_present_label", dtype=ln.ULabel).save()
    missing = ln.Feature(
        name="pl_added_label",
        dtype=ln.ULabel,
        nullable=True,
        default_value="unregistered added label",
    ).save()
    schema = ln.Schema(
        features=[present, missing.with_config(optional=True)],
    ).save()
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_present_label": [known.name]})),
        schema,
    )
    curator.validate()
    curator.standardize()
    assert "pl_added_label" in curator.dataset.collect_schema().names()
    with pytest.raises(ValidationError, match="unregistered added label"):
        curator.validate()


def test_polars_csv_artifact(as_frame):
    schema = make_schema([("pl_csv", int, {"coerce": True})])
    frame = as_frame(pl.DataFrame({"pl_csv": ["1", "2"]}))
    artifact = ln.Artifact.from_dataframe(
        frame, schema=schema, key="polars/coerced.csv"
    ).save()
    assert pl.read_csv(artifact.path)["pl_csv"].to_list() == [1, 2]


def test_polars_coercion_persisted(as_frame):
    schema = make_schema([("pl_persisted", int, {"coerce": True})])
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_persisted": ["1", "2"]})), schema
    )
    artifact = curator.save_artifact(key="polars/coerced.parquet")
    assert pl.read_parquet(artifact.path)["pl_persisted"].dtype == pl.Int64
    assert isinstance(curator.dataset, type(as_frame(pl.DataFrame())))


def test_polars_attrs_not_supported(as_frame):
    attrs = ln.Schema(itype=ln.Feature).save()
    schema = make_schema(
        [("pl_no_attrs", int, {})], slots={"attrs": attrs}, otype="DataFrame"
    )
    with pytest.raises(InvalidArgument, match="do not have an attrs slot"):
        ln.curators.DataFrameCurator(
            as_frame(pl.DataFrame({"pl_no_attrs": [1]})), schema
        )


def test_polars_flexible_schema_is_rejected(as_frame):
    ln.Feature(name="pl_flexible", dtype=int).save()
    schema = ln.Schema(itype=ln.Feature).save()
    with pytest.raises(InvalidArgument, match="flexible schema"):
        ln.curators.DataFrameCurator(
            as_frame(pl.DataFrame({"pl_flexible": [1], "unregistered_column": ["a"]})),
            schema,
        )


def test_polars_schema_coerce_with_check_columns(as_frame):
    schema = make_schema(
        [("pl_cast_int", int, {}), ("pl_no_cast_str", str, {})], coerce=True
    )
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_cast_int": ["1"], "pl_no_cast_str": ["a"]})),
        schema,
    )
    curator.validate()
    assert curator.dataset.collect_schema()["pl_cast_int"] == pl.Int64
    assert curator.dataset.collect_schema()["pl_no_cast_str"] == pl.String


def test_polars_column_rules_and_parquet_sorting_metadata(as_frame):
    import pyarrow.parquet as pq

    schema = make_schema(
        [
            ("pl_rule_id", str, {"unique": True}),
            ("pl_rule_group", str, {"nullable": True}),
            ("pl_rule_rank", int, {}),
        ],
        unique=["pl_rule_group", "pl_rule_rank"],
        sorted_by=["pl_rule_group", ("pl_rule_rank", False)],
    )
    df = pl.DataFrame(
        {
            "pl_rule_id": [None, "ID-2", "ID-3", "ID-4"],
            "pl_rule_group": ["a", "a", "b", "b"],
            "pl_rule_rank": [2, 1, 2, 1],
        }
    )
    curator = ln.curators.DataFrameCurator(as_frame(df), schema)
    curator.validate()
    assert isinstance(curator.dataset, type(as_frame(df)))
    artifact = curator.save_artifact(
        key=f"polars/rules-{type(curator.dataset).__name__}.parquet",
        write_parquet_sorting_metadata=True,
    )
    metadata = pq.ParquetFile(artifact.path).metadata
    assert list(metadata.row_group(0).sorting_columns) == [
        pq.SortingColumn(1, descending=False, nulls_first=False),
        pq.SortingColumn(2, descending=True, nulls_first=False),
    ]
    artifact.delete(permanent=True)

    invalid_cases = [
        (
            df.with_columns(pl.Series("pl_rule_id", ["ID-1", "ID-1", None, "ID-4"])),
            r"\['pl_rule_id'\] must be unique",
        ),
        (
            df.with_columns(pl.Series("pl_rule_id", [None, "ID-2", None, "ID-4"])),
            r"\['pl_rule_id'\] must be unique",
        ),
        (
            df.with_columns(pl.Series("pl_rule_rank", [2, 2, 2, 1])),
            r"\['pl_rule_group', 'pl_rule_rank'\] must be unique",
        ),
        (
            df.with_columns(
                pl.Series("pl_rule_group", ["a", "a", None, None]),
                pl.Series("pl_rule_rank", [2, 1, 2, 2]),
            ),
            r"\['pl_rule_group', 'pl_rule_rank'\] must be unique",
        ),
        (df.with_columns(pl.Series("pl_rule_rank", [1, 2, 2, 1])), "not sorted"),
    ]
    for invalid, match in invalid_cases:
        with pytest.raises(ValidationError, match=match):
            ln.curators.DataFrameCurator(as_frame(invalid), schema).validate()


@pytest.mark.parametrize("backend", ["pandas", "polars", "lazyframe"])
@pytest.mark.parametrize(
    "rows, expected",
    [
        ([], True),
        ([(1.0, 1)], True),
        ([(1.0, 2), (1.0, 2), (1.0, 1), (2.0, 3), (None, 2), (float("nan"), 1)], True),
        ([(1.0, 1), (1.0, 2)], False),
        ([(2.0, 3), (1.0, 2)], False),
        ([(None, 2), (1.0, 1)], False),
        ([(None, 1), (float("nan"), 2)], False),
    ],
)
def test_sort_order_adjacent_rows(backend, rows, expected):
    import pandas as pd
    from lamindb.curators.core import _is_sorted

    if backend == "pandas":
        frame = pd.DataFrame(rows, columns=["a", "b"], index=[0] * len(rows))
    else:
        frame = pl.DataFrame(
            rows, schema={"a": pl.Float64, "b": pl.Int64}, orient="row"
        )
        if backend == "lazyframe":
            frame = frame.lazy()
    assert (
        _is_sorted(frame, [("a", True), ("b", False)], backend != "pandas") is expected
    )


def test_lazyframe_default_sorted_save_uses_native_sink(monkeypatch):
    import pyarrow.parquet as pq

    schema = make_schema([("lazy_sort_value", int, {})], sorted_by=["lazy_sort_value"])
    curator = ln.curators.DataFrameCurator(
        pl.DataFrame({"lazy_sort_value": [1, 2, 3]}).lazy(), schema
    )
    calls = []
    sink = pl.LazyFrame.sink_parquet

    def record_sink(self, path, **kwargs):
        calls.append(path)
        return sink(self, path, **kwargs)

    def no_rewrite(*args, **kwargs):
        raise AssertionError("Default save must not rewrite the Parquet file")

    monkeypatch.setattr(pl.LazyFrame, "sink_parquet", record_sink)
    monkeypatch.setattr(pq, "ParquetWriter", no_rewrite)
    artifact = curator.save_artifact(key="polars/native-sorted.parquet")
    assert len(calls) == 1
    assert not pq.ParquetFile(artifact.path).metadata.row_group(0).sorting_columns
    artifact.delete(permanent=True)


@pytest.mark.parametrize("descending", [False, True])
def test_polars_sort_metadata_nulls_last(as_frame, descending):
    import pyarrow.parquet as pq

    schema = make_schema(
        [("pl_nullable_sort", int, {"nullable": True})],
        sorted_by=[("pl_nullable_sort", not descending)],
    )
    values = [2, 1, None] if descending else [1, 2, None]
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_nullable_sort": values})), schema
    )
    artifact = curator.save_artifact(
        key=f"polars/null-sort-{descending}-{type(curator.dataset).__name__}.parquet",
        write_parquet_sorting_metadata=True,
    )
    metadata = pq.ParquetFile(artifact.path).metadata
    assert list(metadata.row_group(0).sorting_columns) == [
        pq.SortingColumn(0, descending=descending, nulls_first=False)
    ]
    artifact.delete(permanent=True)
    with pytest.raises(ValidationError, match="not sorted"):
        ln.curators.DataFrameCurator(
            as_frame(pl.DataFrame({"pl_nullable_sort": [None, 1, 2]})), schema
        ).validate()


def test_lazyframe_metadata_opt_in_multiple_row_groups():
    import pyarrow.parquet as pq

    schema = make_schema(
        [("lazy_row_group_sort", int, {"nullable": True})],
        sorted_by=["lazy_row_group_sort"],
    )
    df = pl.DataFrame({"lazy_row_group_sort": [1, 2, 3, None]}).lazy()
    artifact = ln.Artifact.from_dataframe(
        df,
        key="polars/multiple-sort-row-groups.parquet",
        schema=schema,
        parquet_kwargs={"row_group_size": 2},
        write_parquet_sorting_metadata=True,
    ).save()
    metadata = pq.ParquetFile(artifact.path).metadata
    assert metadata.num_row_groups == 2
    for i in range(metadata.num_row_groups):
        assert list(metadata.row_group(i).sorting_columns) == [
            pq.SortingColumn(0, descending=False, nulls_first=False)
        ]
    artifact.delete(permanent=True)


def test_lazyframe_metadata_opt_in_streams_batches(monkeypatch, tmp_path):
    import pyarrow.parquet as pq
    from lamindb.core.storage.objects import write_to_disk

    batches = pl.LazyFrame.collect_batches

    def unexpected_collect(*args, **kwargs):
        raise AssertionError("Sorted Parquet writes must not collect the LazyFrame")

    monkeypatch.setattr(
        pl.LazyFrame,
        "collect_batches",
        lambda self, **kwargs: batches(self, chunk_size=2),
    )
    monkeypatch.setattr(pl.LazyFrame, "collect", unexpected_collect)
    sorted_path = tmp_path / "sorted.parquet"
    write_to_disk(
        pl.LazyFrame({"value": [1, 2, 2, 3, None]}),
        sorted_path,
        _lamindb_sorting_columns=[("value", True)],
    )
    metadata = pq.ParquetFile(sorted_path).metadata
    assert metadata.num_row_groups == 3
    assert list(metadata.row_group(2).sorting_columns) == [
        pq.SortingColumn(0, descending=False, nulls_first=False)
    ]
    assert pq.read_table(sorted_path).column("value").to_pylist() == [
        1,
        2,
        2,
        3,
        None,
    ]
    # Each batch is sorted, but the order breaks between batches.
    unsorted_path = tmp_path / "unsorted.parquet"
    with pytest.raises(ValidationError, match="not sorted"):
        write_to_disk(
            pl.LazyFrame({"value": [1, 3, 2, 4]}),
            unsorted_path,
            _lamindb_sorting_columns=[("value", True)],
        )
    assert not unsorted_path.exists()


@pytest.mark.parametrize(
    "dtype,values,lower,upper",
    [
        (int, [0, 5, 10, None], 0, 10),
        (float, [0.0, 5.5, 10.0, float("nan")], 0.0, 10.0),
        ("num", [0.0, 5.5, 10.0, None], 0, 10),
    ],
)
def test_polars_numeric_feature_bounds(as_frame, dtype, values, lower, upper):
    schema = make_schema(
        [
            (
                "pl_bounded_numeric",
                dtype,
                {"min_value": lower, "max_value": upper},
            )
        ]
    )
    ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_bounded_numeric": values})), schema
    ).validate()
    for invalid_values, match in (
        ([lower - 1, upper], "min_value"),
        ([lower, upper + 1], "max_value"),
    ):
        with pytest.raises(ValidationError, match=match):
            ln.curators.DataFrameCurator(
                as_frame(pl.DataFrame({"pl_bounded_numeric": invalid_values})),
                schema,
            ).validate()


def test_polars_external_features(as_frame):
    external = ln.Feature(name="pl_external", dtype=str).save()
    external_schema = ln.Schema(features=[external]).save()
    schema = make_schema(
        [("pl_internal", int, {})],
        slots={"__external__": external_schema},
        otype="DataFrame",
    )
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_internal": [1]})),
        schema,
        features={"pl_external": "external value"},
    )
    artifact = curator.save_artifact(key="polars/external.parquet")
    assert artifact.features.get_values()["pl_external"] == "external value"


def test_polars_null_lists_do_not_match_scalar_string(as_frame):
    schema = make_schema([("pl_wrong_list", str, {})])
    with pytest.raises(ValidationError):
        ln.curators.DataFrameCurator(
            as_frame(pl.DataFrame({"pl_wrong_list": [[], []]})), schema
        ).validate()


def test_polars_standardize_empty_frame_preserves_rows(as_frame):
    schema = make_schema([("pl_empty_default", str, {"default_value": "filled"})])
    curator = ln.curators.DataFrameCurator(as_frame(pl.DataFrame()), schema)
    curator.standardize()
    assert eager(curator.dataset).height == 0
    curator.validate()


@pytest.mark.parametrize("suffix", [".parquet", ".csv", ".tsv"])
def test_from_dataframe_path_validates_lazily(monkeypatch, tmp_path, suffix):
    schema = make_schema([("pl_file", int, {})])
    path = tmp_path / f"file{suffix}"
    frame = pl.DataFrame({"pl_file": [1, 2, 3]})
    if suffix == ".parquet":
        frame.write_parquet(path)
    else:
        frame.write_csv(path, separator="\t" if suffix == ".tsv" else ",")

    def unexpected_load(*args, **kwargs):
        pytest.fail("The file must be scanned lazily, not loaded")

    with monkeypatch.context() as context:
        context.setattr(pl.LazyFrame, "collect", unexpected_load)
        context.setattr(ln.Artifact, "load", unexpected_load)
        artifact = ln.Artifact.from_dataframe(
            str(path), key=f"polars/scan{suffix}", schema=schema
        )
    assert isinstance(artifact._curator.dataset, pl.LazyFrame)
    artifact.save()
    assert artifact.schema == schema
    assert artifact.path.read_bytes() == path.read_bytes()

    bad = tmp_path / f"bad{suffix}"
    bad_frame = pl.DataFrame({"pl_file": ["x"]})
    if suffix == ".parquet":
        bad_frame.write_parquet(bad)
    else:
        bad_frame.write_csv(bad, separator="\t" if suffix == ".tsv" else ",")
    with pytest.raises(ValidationError):
        ln.Artifact.from_dataframe(str(bad), key=f"polars/bad{suffix}", schema=schema)


def test_registry_checks_use_one_batched_query(as_frame):
    labels = [ln.ULabel(name=f"pl_batch_{i}").save() for i in range(2)]
    schema = make_schema(
        [
            ("pl_batch_a", ln.ULabel, {}),
            ("pl_batch_b", ln.ULabel, {}),
        ]
    )
    names = [label.name for label in labels]
    frame = as_frame(
        pl.DataFrame(
            {
                "pl_batch_a": names * 50,
                "pl_batch_b": names[::-1] * 50,
            }
        )
    )
    curator = ln.curators.DataFrameCurator(frame, schema)
    curator.validate()
    manager = curator._atomic_curator.cat
    values = manager._get_unique_values()
    assert sorted(values["pl_batch_a"]) == names
    assert len(values["pl_batch_a"]) == 2


@pytest.mark.parametrize("suffix", [".parquet", ".csv"])
def test_from_dataframe_path_stores_coerced_values(tmp_path, suffix):
    schema = make_schema([("pl_file_coerce", int, {"coerce": True})])
    path = tmp_path / f"coerce{suffix}"
    frame = pl.DataFrame({"pl_file_coerce": ["1", "2"]})
    # string dtype is kept in the source file, the CSV reader would infer ints
    frame.write_parquet(path) if suffix == ".parquet" else frame.write_csv(path)
    artifact = ln.Artifact.from_dataframe(
        str(path), key=f"polars/coerce{suffix}", schema=schema
    ).save()
    assert artifact.suffix == suffix
    assert (
        pl.scan_parquet(artifact.path).collect_schema()["pl_file_coerce"] == pl.Int64
        if suffix == ".parquet"
        else True
    )
    stored = (
        pl.read_parquet(artifact.path)
        if suffix == ".parquet"
        else pl.read_csv(artifact.path)
    )
    assert stored["pl_file_coerce"].to_list() == [1, 2]
    assert stored["pl_file_coerce"].dtype == pl.Int64


def test_flexible_schema_path_uses_pandas(tmp_path):
    ln.Feature(name="pl_flexible_file", dtype=int).save()
    schema = ln.Schema(itype=ln.Feature).save()
    path = tmp_path / "flexible.parquet"
    pl.DataFrame({"pl_flexible_file": [1], "extra": ["a"]}).write_parquet(path)
    artifact = ln.Artifact.from_dataframe(
        str(path), key="polars/flexible.parquet", schema=schema
    )
    assert not isinstance(artifact._curator.dataset, (pl.DataFrame, pl.LazyFrame))
    artifact.save()
    assert artifact.schema == schema
    assert artifact.path.read_bytes() == path.read_bytes()


def test_scan_dataframe_file_falls_back_on_failure(tmp_path):
    from lamindb.models.artifact import _scan_dataframe_file

    assert _scan_dataframe_file(UPath(tmp_path / "missing.parquet")) is None


def test_remote_paths_are_scannable(monkeypatch):
    from lamindb.models.artifact import _can_scan_lazily

    schema = make_schema([("pl_remote", int, {})])
    assert _can_scan_lazily("s3://bucket/data.parquet", schema)
    assert _can_scan_lazily("gs://bucket/data.csv", schema)
    assert _can_scan_lazily("s3://bucket/data.tsv", schema)
    assert not _can_scan_lazily("s3://bucket/data.h5ad", schema)
    assert not _can_scan_lazily("ftp://host/data.parquet", schema)

    coercing = make_schema([("pl_remote_coerce", int, {"coerce": True})])
    assert not _can_scan_lazily("s3://bucket/data.tsv", coercing)

    index = ln.Feature(name="pl_scan_index", dtype=str).save()
    indexed = ln.Schema(itype=ln.Feature, index=index).save()
    assert not _can_scan_lazily("s3://bucket/data.parquet", indexed)

    flexible = ln.Schema(itype=ln.Feature).save()
    assert not _can_scan_lazily("s3://bucket/data.parquet", flexible)

    monkeypatch.setattr(
        "lamindb.models.artifact.is_package_installed", lambda _name: False
    )
    assert not _can_scan_lazily("s3://bucket/data.parquet", schema)
