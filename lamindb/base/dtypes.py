"""Dtype utils.

.. autofunction:: check_dtype

"""

from datetime import datetime
from typing import Any, Callable, Iterable

import pandas as pd


def is_list_of_type(value: Any, expected_type: Any) -> bool:
    """Helper function to check if a value is either of expected_type or a list of that type, or a mix of both in a nested structure."""
    if isinstance(value, Iterable) and not isinstance(value, (str, bytes)):
        # handle nested lists recursively
        return all(isinstance(item, expected_type) for item in value)
    return False


def check_dtype(expected_type: Any, nullable: bool) -> Callable:
    """Creates a check function for Pandera that validates a column's dtype.

    Supports both standard dtype checking and mixed list/single values for the same type.
    For example, a column with expected_type 'float' would also accept a mix of float values and lists of floats.

    Args:
        expected_type: String identifier for the expected type ('int', 'float', 'num', 'str')

    Returns:
        A function that checks if a series has the expected dtype or contains mixed types
    """

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
        elif expected_type == "str" and pd.api.types.is_string_dtype(series.dtype):
            return True
        elif expected_type == "path" and pd.api.types.is_string_dtype(series.dtype):
            return True

        # if we're here, it might be a mixed column with object dtype
        # need to check each value individually
        if series.dtype == "object" and expected_type.startswith("list"):
            expected_type_member = expected_type.replace("list[", "").removesuffix("]")
            if expected_type_member == "int":
                return series.apply(lambda x: is_list_of_type(x, int)).all()
            elif expected_type_member == "float":
                return series.apply(lambda x: is_list_of_type(x, float)).all()
            elif expected_type_member == "num":
                # for numeric, accept either int or float
                return series.apply(lambda x: is_list_of_type(x, (int, float))).all()
            elif (
                expected_type_member == "str"
                or expected_type_member == "path"
                or expected_type_member.startswith("cat[")
            ):
                return series.apply(lambda x: is_list_of_type(x, str)).all()
            elif expected_type_member == "list":
                return series.apply(lambda x: isinstance(x, list)).all()

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
