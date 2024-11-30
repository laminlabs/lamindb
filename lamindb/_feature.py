from __future__ import annotations

from typing import TYPE_CHECKING, Any, Literal, get_args

import lamindb_setup as ln_setup
import pandas as pd
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from lnschema_core.models import Artifact, Feature, Record
from lnschema_core.types import FeatureDtype
from pandas.api.types import CategoricalDtype, is_string_dtype

from lamindb.core.exceptions import ValidationError

from ._query_set import RecordList
from ._utils import attach_func_to_class_method
from .core._settings import settings
from .core.schema import dict_schema_name_to_model_name

if TYPE_CHECKING:
    from collections.abc import Iterable

    from lnschema_core.types import FieldAttr
    from pandas.core.dtypes.base import ExtensionDtype


FEATURE_DTYPES = set(get_args(FeatureDtype))


def get_dtype_str_from_dtype(dtype: Any) -> str:
    if not isinstance(dtype, list) and dtype.__name__ in FEATURE_DTYPES:
        dtype_str = dtype.__name__
    else:
        error_message = "dtype has to be of type Record or list[Record]"
        if isinstance(dtype, Record):
            dtype = [dtype]
        elif not isinstance(dtype, list):
            raise ValueError(error_message)
        registries_str = ""
        for registry in dtype:
            if not hasattr(registry, "__get_name_with_schema__"):
                raise ValueError(error_message)
            registries_str += registry.__get_name_with_schema__() + "|"
        dtype_str = f'cat[{registries_str.rstrip("|")}]'
    return dtype_str


def convert_pandas_dtype_to_lamin_dtype(pandas_dtype: ExtensionDtype) -> str:
    if is_string_dtype(pandas_dtype):
        if not isinstance(pandas_dtype, CategoricalDtype):
            dtype = "str"
        else:
            dtype = "cat"
    # there are string-like categoricals and "pure" categoricals (pd.Categorical)
    elif isinstance(pandas_dtype, CategoricalDtype):
        dtype = "cat"
    else:
        # strip precision qualifiers
        dtype = "".join(dt for dt in pandas_dtype.name if not dt.isdigit())
    if dtype.startswith("datetime"):
        dtype = dtype.split("[")[0]
    assert dtype in FEATURE_DTYPES  # noqa: S101
    return dtype


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(Feature, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) != 0:
        raise ValueError("Only keyword args allowed")
    dtype: type | str = kwargs.pop("dtype") if "dtype" in kwargs else None
    # cast type
    if dtype is None:
        raise ValueError(f"Please pass dtype, one of {FEATURE_DTYPES}")
    elif dtype is not None:
        if not isinstance(dtype, str):
            dtype_str = get_dtype_str_from_dtype(dtype)
        else:
            dtype_str = dtype
            # add validation that a registry actually exists
            if dtype_str not in FEATURE_DTYPES and not dtype_str.startswith("cat"):
                raise ValueError(
                    f"dtype is {dtype_str} but has to be one of {FEATURE_DTYPES}!"
                )
            if dtype_str != "cat" and dtype_str.startswith("cat"):
                registries_str = dtype_str.replace("cat[", "").rstrip("]")
                if registries_str != "":
                    registry_str_list = registries_str.split("|")
                    for registry_str in registry_str_list:
                        if registry_str not in dict_schema_name_to_model_name(Artifact):
                            raise ValueError(
                                f"'{registry_str}' is an invalid dtype, pass, e.g. `[ln.ULabel, bt.CellType]` or similar"
                            )
    kwargs["dtype"] = dtype_str
    super(Feature, self).__init__(*args, **kwargs)
    if not self._state.adding:
        if not (
            self.dtype.startswith("cat") if dtype == "cat" else self.dtype == dtype
        ):
            raise ValidationError(
                f"Feature {self.name} already exists with dtype {self.dtype}, you passed {dtype}"
            )


def suggest_categorical_for_str_iterable(
    iterable: Iterable[str], key: str = None
) -> str:
    c = pd.Categorical(iterable)
    message = ""
    if len(c.categories) < len(c):
        if key != "":
            key_note = f" for feature {key}"
        else:
            key_note = ""
        message = f"You have few permissible values{key_note}, consider dtype 'cat' instead of 'str'"
    return message


def categoricals_from_df(df: pd.DataFrame) -> dict:
    """Returns categorical columns."""
    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {
        col: df[col]
        for col in df.columns
        if isinstance(df[col].dtype, CategoricalDtype)
    }
    for key in string_cols:
        message = suggest_categorical_for_str_iterable(df[key], key)
        if message:
            logger.warning(message)
    return categoricals


@classmethod  # type:ignore
@doc_args(Feature.from_df.__doc__)
def from_df(cls, df: pd.DataFrame, field: FieldAttr | None = None) -> RecordList:
    """{}"""  # noqa: D415
    field = Feature.name if field is None else field
    registry = field.field.model
    if registry != Feature:
        raise ValueError("field must be a Feature FieldAttr!")
    categoricals = categoricals_from_df(df)
    dtypes = {}
    for name, col in df.items():
        if name in categoricals:
            dtypes[name] = "cat"
        else:
            dtypes[name] = convert_pandas_dtype_to_lamin_dtype(col.dtype)
    with logger.mute():  # silence the warning "loaded record with exact same name "
        features = [Feature(name=name, dtype=dtype) for name, dtype in dtypes.items()]
    assert len(features) == len(df.columns)  # noqa: S101
    return RecordList(features)


@doc_args(Feature.save.__doc__)
def save(self, *args, **kwargs) -> Feature:
    """{}"""  # noqa: D415
    super(Feature, self).save(*args, **kwargs)
    return self


METHOD_NAMES = [
    "__init__",
    "from_df",
    "save",
]

if ln_setup._TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(Feature, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Feature, globals())
