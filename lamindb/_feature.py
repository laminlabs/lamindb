from __future__ import annotations

from typing import TYPE_CHECKING, Any

import lamindb_setup as ln_setup
import pandas as pd
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from lnschema_core.models import Artifact, Feature, Record
from pandas.api.types import CategoricalDtype, is_string_dtype

from lamindb.core.exceptions import ValidationError

from ._query_set import RecordsList
from ._utils import attach_func_to_class_method
from .core._settings import settings
from .core.schema import dict_schema_name_to_model_name

if TYPE_CHECKING:
    from lnschema_core.types import FieldAttr


FEATURE_DTYPES = {
    "cat",  # categorical variables
    "num",  # numerical variables
    "str",  # string variables
    "int",  # integer variables
    "float",  # float variables
    "bool",  # boolean variables
    "date",  # date variables
    "datetime",  # datetime variables
    "object",  # this is a pandas type, we're only using it for complicated types, not for strings
}


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


def convert_pandas_dtype_to_lamin_dtype(pandas_dtype) -> str:
    if is_string_dtype(pandas_dtype) and not isinstance(pandas_dtype, CategoricalDtype):
        dtype = "str"
    else:
        # strip precision qualifiers
        dtype = "".join(i for i in pandas_dtype.name if not i.isdigit())
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
        raise ValueError("Please pass dtype!")
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
        if self.dtype != dtype:
            raise ValidationError(
                f"Feature already exists with dtype {self.dtype}, you passed {dtype}"
            )


def categoricals_from_df(df: pd.DataFrame) -> dict:
    """Returns categorical columns."""
    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {
        col: df[col]
        for col in df.columns
        if isinstance(df[col].dtype, CategoricalDtype)
    }
    for key in string_cols:
        c = pd.Categorical(df[key])
        if len(c.categories) < len(c):
            logger.warning(
                "consider changing the dtype of string column `key` to categorical"
            )
    return categoricals


@classmethod  # type:ignore
@doc_args(Feature.from_df.__doc__)
def from_df(cls, df: pd.DataFrame, field: FieldAttr | None = None) -> RecordsList:
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
        # create records for all features including non-validated
        features = [Feature(name=name, dtype=dtype) for name, dtype in dtypes.items()]
    assert len(features) == len(df.columns)  # noqa: S101
    return RecordsList(features)


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
