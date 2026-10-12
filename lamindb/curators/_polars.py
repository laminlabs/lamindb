"""Native Pandera validation for the optional Polars backend."""

import polars as pl
from pandera import polars as pandera
from pandera.constants import CHECK_OUTPUT_KEY
from pandera.engines import polars_engine
from pandera.errors import ParserError


def _all_missing(data) -> bool:
    return data.lazyframe.select(pl.col(data.key).is_null().all()).collect().item()


class _NumericCoercion:
    def try_coerce(self, data_container):
        coerced = self.coerce(data_container)
        try:
            # Evaluate the cast, not the entire LazyFrame.
            coerced.select(pl.col(data_container.key).is_null().sum()).collect()
        except (pl.exceptions.InvalidOperationError, pl.exceptions.ComputeError) as err:
            values = pl.col(data_container.key)
            valid = (
                values.is_null() | values.cast(self.type, strict=False).is_not_null()
            )
            raise ParserError(
                f"Could not coerce into {self.type}",
                failure_cases=data_container.lazyframe.filter(~valid)
                .select(data_container.key)
                .collect(),
                parser_output=data_container.lazyframe.select(
                    valid.alias(CHECK_OUTPUT_KEY)
                ).collect(),
            ) from err
        return coerced


@polars_engine.Engine.register_dtype(equivalents=["lamindb.int"])
@polars_engine.immutable
class AnyInt(_NumericCoercion, polars_engine.Int64):
    """Accept all integer widths and prevent truncating float coercions."""

    def check(self, pandera_dtype, data_container=None):
        dtype = polars_engine.Engine.dtype(pandera_dtype).type
        return dtype.is_integer() or (
            data_container is not None and _all_missing(data_container)
        )

    def coerce(self, data_container):
        dtype = data_container.lazyframe.collect_schema()[data_container.key]
        if dtype.is_integer():
            return data_container.lazyframe
        if dtype.is_float():
            values = pl.col(data_container.key)
            valid = values.is_null() | (values == values.floor())
            lossless = data_container.lazyframe.select(valid.all()).collect().item()
            if not lossless:
                raise ParserError(
                    "Could not losslessly coerce into int",
                    failure_cases=data_container.lazyframe.filter(~valid)
                    .select(data_container.key)
                    .collect(),
                    parser_output=data_container.lazyframe.select(
                        valid.alias(CHECK_OUTPUT_KEY)
                    ).collect(),
                )
        return super().coerce(data_container)


@polars_engine.Engine.register_dtype(equivalents=["lamindb.float"])
@polars_engine.immutable
class AnyFloat(_NumericCoercion, polars_engine.Float64):
    """Accept all float widths, preserving existing floating-point dtypes."""

    def check(self, pandera_dtype, data_container=None):
        dtype = polars_engine.Engine.dtype(pandera_dtype).type
        return dtype.is_float() or (
            data_container is not None and _all_missing(data_container)
        )

    def coerce(self, data_container):
        dtype = data_container.lazyframe.collect_schema()[data_container.key]
        if dtype.is_float():
            return data_container.lazyframe
        return super().coerce(data_container)


def _matches_dtype(dtype, expected: str) -> bool:
    if expected.startswith("list["):
        return isinstance(dtype, pl.List) and _matches_dtype(
            dtype.inner, expected[5:-1]
        )
    if expected == "list":
        return isinstance(dtype, pl.List)
    if expected in {"str", "path", "url"} or expected.startswith("cat"):
        return dtype == pl.String or isinstance(dtype, (pl.Categorical, pl.Enum))
    if expected == "int":
        return dtype.is_integer()
    if expected == "float":
        return dtype.is_float()
    if expected == "num":
        return dtype.is_integer() or dtype.is_float()
    if expected == "bool":
        return dtype == pl.Boolean
    if expected == "dict":
        return isinstance(dtype, pl.Struct)
    return False


def column(feature, required: bool, schema_coerce: bool = False):
    """Translate a Lamin feature into a native Polars Pandera column."""
    dtype_str = feature._dtype_str
    kwargs = {
        "nullable": feature.nullable,
        "coerce": feature.coerce or schema_coerce,
        "required": required,
    }
    if dtype_str in {"int", "float"}:
        return pandera.Column(f"lamindb.{dtype_str}", **kwargs)
    if (
        dtype_str in {"bool", "num", "str", "path", "url", "dict"}
        or dtype_str.startswith("list")
        or dtype_str.startswith("cat")
    ):
        kwargs["coerce"] = False

        def check(data):
            dtype = data.lazyframe.collect_schema()[data.key]
            if (
                dtype_str.startswith("list")
                and isinstance(dtype, pl.List)
                and dtype.inner == pl.Null
            ):
                return (
                    data.lazyframe.select(
                        (
                            pl.col(data.key).is_null()
                            | (pl.col(data.key).list.len() == 0)
                        ).all()
                    )
                    .collect()
                    .item()
                )
            if dtype == pl.Null and (
                dtype_str == "str" or dtype_str.startswith("list")
            ):
                return True
            return _matches_dtype(dtype, dtype_str) or (
                feature.nullable and _all_missing(data)
            )

        return pandera.Column(
            dtype=None,
            checks=pandera.Check(
                check,
                error=f"Column '{feature.name}' failed dtype check for '{dtype_str}'",
            ),
            **kwargs,
        )
    dtype = (
        pl.Datetime(time_unit="ns", time_zone="UTC")
        if dtype_str == "datetime64[ns, UTC]"
        else pl.Datetime
        if dtype_str == "datetime"
        else pl.Date
        if dtype_str == "date"
        else dtype_str
    )
    return pandera.Column(dtype, **kwargs)
