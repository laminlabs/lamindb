"""Data types.

Most data types in LaminDB are string-serializable representations of standard Python types, commonly used for validating data in libraries such as pydantic.
Categorical data types are absent from built-in Python types and are needed to anchor the validation of a categorical variable in a LaminDB registry.

Simple
------

In the table below, the first column shows the object that
can be passed to the `dtype` argument of `Feature()` or `Schema()` and the second the string serialization
that's used in the database.

.. list-table::
    :header-rows: 1

    * - dtype
      - string serialization
      - pandas
    * - `int`
      - `"int"`
      - `int64 | int32 | int16 | int8 | uint | ...`
    * - `float`
      - `"float"`
      - `float64 | float32 | float16 | float8 | ...`
    * - `str`
      - `"str"`
      - `object`
    * - `bool`
      - `"bool"`
      - `boolean | bool`
    * - `datetime`
      - `"datetime"`
      - `datetime`
    * - `"datetime64[ns, UTC]"`
      - `"datetime64[ns, UTC]"`
      - `datetime64[ns, UTC]`
    * - `date`
      - `"date"`
      - `object` (pandera requires an ISO-format string, convert with `df["date"] = df["date"].dt.date`)
    * - `dict`
      - `"dict"`
      - `object`
    * - `"num"`
      - `"num"`
      - `int | float` ("num" is a convenience type for `int | float`)
    * - `"path"`
      - `"path"`
      - `str` (pandas does not have a dedicated path type, validated as `str`)
    * - `"url"`
      - `"url"`
      - `str` (pandas does not have a dedicated url type, validated as `str`)

Categorical/ relational
-----------------------

For any categorical, you can restrict permissible
values to the values defined in a registry. This establishes a relationship.

.. list-table::
    :header-rows: 1

    * - dtype
      - string serialization
    * - `ln.ULabel`
      - `"cat[ULabel]"`
    * - `bt.CellType`
      - `"cat[bionty.CellType]"`
    * - `bt.Disease`
      - `"cat[bionty.Disease]"`
    * - `ln.Artifact`
      - `"cat[Artifact]"`

You can restrict permissible values to instances of `ULabel` or `Record` types, i.e., to dynamic registries.

.. list-table::
    :header-rows: 1

    * - dtype
      - string serialization
    * - `ulabel_type` (a `ULabel` with `is_type=True`)
      - `"cat[ULabel[<uid_of_ulabel_type>]]"`
    * - `record_type` (a `Record` with `is_type=True`)
      - `"cat[Record[<uid_of_record_type>]]"`

You can restrict permissible values by filtering the categorical on fields of its registry.

.. list-table::
    :header-rows: 1

    * - dtype
      - cat_filters
      - string serialization
    * - `bt.Disease`
      - `{"source": source}`
      - `"cat[bionty.Disease[source__uid='<uid_of_source>']]"`
    * - `ln.Artifact`
      - `{"schema": schema}`
      - `"cat[Artifact[schema__uid='<uid_of_schema>']]"`

Lists
-----

.. list-table::
    :header-rows: 1

    * - dtype
      - string serialization
    * - `list[bt.CellType]`
      - `"list[cat[bionty.CellType]]"`
    * - `list[float]`
      - `"list[float]"`

Unions
------

Unions are currently only supported for static registries.

.. list-table::
    :header-rows: 1

    * - dtype
      - string serialization
    * - `[bt.Tissue.ontology_id, bt.CellType.ontology_id]`
      - `"cat[bionty.Tissue.ontology_id|bionty.CellType.ontology_id]"`

Usage
-----

In features, you can pass dtypes to the `dtype` argument. See :class:`~lamindb.Feature` for more details.

In function signatures, you can use dtypes to type annotate and validate arguments.
See :doc:`docs:track` or :func:`~lamindb.flow` for examples.

"""

from collections.abc import Iterable, Sequence
from datetime import datetime
from typing import Any, Callable

import numpy as np
import pandas as pd
from pandera.engines import pandas_engine


def is_list_of_type(value: Any, expected_type: Any) -> bool:
    """Helper function to check if a value is either of expected_type or a list of that type, or a mix of both in a nested structure."""
    if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        # handle nested lists recursively
        return all(isinstance(item, expected_type) for item in value)
    return False


def check_pandera_str(series) -> bool:
    """Validate a series/index as lamin `str`, matching pandas 2 pandera results.

    Pandera maps `Column("str")` differently by pandas version:

    - pandas 2: `NpString` — default strings are `object`
    - pandas 3: `STRING` / `string[pyarrow]` — default strings are `str`

    Target results (pandas 2 `Column("str")`):

    - `object` all-str → accept
    - `object` mixed → reject (elementwise)
    - `object` empty → accept
    - `string` / `string[pyarrow]` → accept
    - string `category` → accept
    - `int64` / other → reject

    On pandas 3, bare `Engine.dtype("str").check` already matches all-str /
    mixed `object` and string dtypes, but rejects empty `object` and
    `category`. Those two cases are handled explicitly below.
    """
    import pandas as pd
    from pandera.engines import pandas_engine

    # any empty series: pandas 2 NpString accepts vacuously (incl. empty int64);
    # pandas 3 STRING rejects empty object/int. also covers empty export RangeIndex
    if len(series) == 0:
        return True
    # string category: accept on pandas 2 (elementwise), reject on pandas 3 STRING.
    # AnnData often stores str obs as categorical
    if isinstance(series.dtype, pd.CategoricalDtype):
        return all(isinstance(x, str) for x in series.dtype.categories)

    # pandas 2 and 3 agree here: string dtypes accept; object is checked
    # elementwise (all-str accept, mixed reject)
    result = pandas_engine.Engine.dtype("str").check(series.dtype, series)
    # object dtype → iterable of bools; string dtypes → scalar bool
    if isinstance(result, bool):
        return result
    return bool(all(result))


def try_coerce_simple_dtype(series, expected_type: str):
    """Losslessly coerce a Series to `int` or `float`, or return `None`.

    Used by :class:`AnyInt` and :class:`AnyFloat`. Does not truncate
    (e.g. `1.1` → `int` fails). Returns the original series if it already
    has the expected pandas dtype, without changing its width.
    """
    import pandas as pd

    if expected_type == "int":
        if pd.api.types.is_integer_dtype(series.dtype):
            return series
        try:
            numeric = pd.to_numeric(series, errors="raise")
        except (TypeError, ValueError):
            return None
        non_null = numeric.dropna()
        if len(non_null) and not bool((non_null == non_null.round()).all()):
            return None
        if numeric.hasnans:
            return numeric.astype("Int64")
        return numeric.astype("int64")
    if expected_type == "float":
        if pd.api.types.is_float_dtype(series.dtype):
            return series
        try:
            return pd.to_numeric(series, errors="raise").astype("float64")
        except (TypeError, ValueError):
            return None
    return None


def _accepts_any_width(series, kind_check) -> bool:
    """True when every value is missing, or the series dtype passes `kind_check`."""
    if series is None:
        return False
    if len(series) == 0 or bool(series.isna().all()):
        return True
    return bool(kind_check(series.dtype))


def _coerce_lossless(series, expected_type: str):
    """Return a lossless coercion, or raise pandera's ParserError."""
    from pandera.errors import ParserError

    coerced = try_coerce_simple_dtype(series, expected_type)
    if coerced is None:
        raise ParserError(
            f"Could not losslessly coerce into {expected_type}",
            failure_cases=series,
        )
    return coerced


# Registered with no "int"/"float" equivalents, so pandera's fixed-width
# dtypes stay unchanged. Same hook as pandas_engine.DateTime: a real dtype,
# which is what pandera coerces. A Check cannot honor coerce.
@pandas_engine.Engine.register_dtype
@pandas_engine.immutable
class AnyInt(pandas_engine.DataType):
    """Integer dtype that accepts any width and coerces only losslessly."""

    type = np.dtype("int64")

    def check(self, pandera_dtype, data_container=None):
        return _accepts_any_width(data_container, pd.api.types.is_integer_dtype)

    def coerce(self, data_container):
        return _coerce_lossless(data_container, "int")

    def __str__(self) -> str:
        return "int"


@pandas_engine.Engine.register_dtype
@pandas_engine.immutable
class AnyFloat(pandas_engine.DataType):
    """Float dtype that accepts any width and coerces only losslessly."""

    type = np.dtype("float64")

    def check(self, pandera_dtype, data_container=None):
        return _accepts_any_width(data_container, pd.api.types.is_float_dtype)

    def coerce(self, data_container):
        return _coerce_lossless(data_container, "float")

    def __str__(self) -> str:
        return "float"


def check_dtype(expected_type: Any, nullable: bool) -> Callable:
    """Creates a check function for Pandera that validates a column's dtype.

    Supports both standard dtype checking and mixed list/single values for the same type.
    For example, a column with expected_type 'float' would also accept a mix of float values and lists of floats.

    Args:
        expected_type: String identifier for the expected type ('int', 'float', 'num', 'str')

    Returns:
        A function that checks if a series has the expected dtype or contains mixed types
    """
    import pandas as pd

    def check_function(series):
        # empty series are considered valid if feature is nullable
        # the issue is that nullable in Pandera controls whether None/NaN values are allowed in the column, not whether the column can be empty (0 rows).
        # so "col": [1, 2, None, 4] is correctly handled by pandera nullable=True, but an empty column "col": [] is not.
        if nullable and series.isnull().all():
            return True
        # first check if the series is entirely of the expected dtype (fast path)
        if expected_type == "int" and pd.api.types.is_integer_dtype(series.dtype):
            return True
        elif expected_type == "float" and pd.api.types.is_float_dtype(series.dtype):
            return True
        elif expected_type == "num" and pd.api.types.is_numeric_dtype(series.dtype):
            return True
        elif expected_type == "str":
            return check_pandera_str(series)
        elif expected_type == "path" and pd.api.types.is_string_dtype(series.dtype):
            return True
        elif expected_type == "url" and pd.api.types.is_string_dtype(series.dtype):
            return True
        elif expected_type == "bool" and pd.api.types.is_bool_dtype(series.dtype):
            return True

        # if we're here, it might be a mixed column with object dtype
        # need to check each value individually
        if series.dtype == "object" and expected_type.startswith("list"):
            expected_type_member = expected_type.replace("list[", "").removesuffix("]")
            if expected_type_member == "int":
                return series.apply(lambda x: is_list_of_type(x, int)).all()
            elif expected_type_member == "float":
                return series.apply(lambda x: is_list_of_type(x, float)).all()
            elif expected_type_member == "bool":
                return series.apply(lambda x: is_list_of_type(x, bool)).all()
            elif expected_type_member == "num":
                # for numeric, accept either int or float
                return series.apply(lambda x: is_list_of_type(x, (int, float))).all()
            elif (
                expected_type_member == "str"
                or expected_type_member == "path"
                or expected_type_member == "url"
                or expected_type_member.startswith("cat[")
            ):
                return series.apply(lambda x: is_list_of_type(x, str)).all()
            elif expected_type_member == "list":
                return series.apply(
                    lambda x: isinstance(x, Sequence)
                    and not isinstance(x, (str, bytes))
                ).all()

        # if we get here, the validation failed
        return False

    return check_function


def is_valid_datetime_str(date_string: str) -> bool | str:
    try:
        dt = datetime.fromisoformat(date_string)
        return dt.isoformat()
    except ValueError:
        return False


def is_iterable_of_sqlrecord(value: Any):
    from lamindb.models import SQLRecord

    return isinstance(value, Iterable) and isinstance(next(iter(value)), SQLRecord)
