from datetime import date, datetime

import bionty as bt
import lamindb as ln
import pytest
from django.db import transaction
from lamindb.errors import InvalidArgument, ValidationError

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


@pytest.mark.parametrize(
    "dtype,values",
    [
        (int, ["a"]),
        (float, [1]),
        (str, [1]),
        (bool, [1]),
        (list[str], [[1]]),
        (dict, ["a"]),
    ],
)
def test_polars_invalid_dtype(as_frame, dtype, values):
    schema = make_schema([("pl_invalid", dtype, {})])
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_invalid": values})), schema
    )
    with pytest.raises(ValidationError):
        curator.validate()
    assert not curator._is_validated


@pytest.mark.parametrize("coerce_on", ["feature", "schema"])
@pytest.mark.parametrize("dtype,expected", [(int, pl.Int64), (float, pl.Float64)])
def test_polars_coercion(as_frame, coerce_on, dtype, expected):
    schema = make_schema(
        [("pl_coerce", dtype, {"coerce": coerce_on == "feature"})],
        coerce=coerce_on == "schema",
    )
    df = pl.DataFrame({"pl_coerce": ["1", "2"]})
    curator = ln.curators.DataFrameCurator(as_frame(df), schema)
    curator.validate()
    assert type(curator.dataset) is type(as_frame(df))
    assert curator.dataset.collect_schema()["pl_coerce"] == expected
    assert eager(curator.dataset)["pl_coerce"].to_list() == [1, 2]


def test_polars_coercion_is_lossless(as_frame):
    schema = make_schema([("pl_lossless", int, {"coerce": True})])
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_lossless": [1.5]})), schema
    )
    with pytest.raises(ValidationError):
        curator.validate()


@pytest.mark.parametrize("nullable", [True, False])
def test_polars_nullability(as_frame, nullable):
    schema = make_schema([("pl_nullable", str, {"nullable": nullable})])
    curator = ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_nullable": ["a", None]})), schema
    )
    if nullable:
        curator.validate()
    else:
        with pytest.raises(ValidationError):
            curator.validate()


def test_polars_column_constraints(as_frame):
    schema = make_schema(
        [("pl_first", int, {}), ("pl_second", str, {})],
        maximal_set=True,
        ordered_set=True,
    )
    for df in [
        pl.DataFrame({"pl_first": [1]}),
        pl.DataFrame({"pl_second": ["a"], "pl_first": [1]}),
        pl.DataFrame({"pl_first": [1], "pl_second": ["a"], "extra": [1]}),
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
    optional = ln.Feature(name="pl_optional", dtype=str).save()
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


def test_polars_flexible_schema(as_frame):
    ln.Feature(name="pl_flexible", dtype=int).save()
    schema = ln.Schema(itype=ln.Feature).save()
    ln.curators.DataFrameCurator(
        as_frame(pl.DataFrame({"pl_flexible": [1], "unregistered_column": ["a"]})),
        schema,
    ).validate()


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


def test_scan_dataframe_file_falls_back_on_failure(tmp_path):
    from lamindb.models.artifact import _scan_dataframe_file

    assert _scan_dataframe_file(str(tmp_path / "missing.parquet")) is None


def test_remote_paths_are_scannable():
    from lamindb.models.artifact import _can_scan_lazily

    schema = make_schema([("pl_remote", int, {})])
    assert _can_scan_lazily("s3://bucket/data.parquet", schema)
    assert _can_scan_lazily("gs://bucket/data.csv", schema)
    assert not _can_scan_lazily("s3://bucket/data.h5ad", schema)
    assert not _can_scan_lazily("ftp://host/data.parquet", schema)
