from __future__ import annotations

from typing import TYPE_CHECKING

import lamindb_setup as ln_setup
import pandas as pd
from lamindb_setup.core._docs import doc_args
from lnschema_core.models import Artifact, Feature
from pandas.api.types import CategoricalDtype, is_string_dtype

from lamindb._utils import attach_func_to_class_method
from lamindb.core._settings import settings

from ._query_set import RecordsList
from .core.schema import dict_schema_name_to_model_name

if TYPE_CHECKING:
    from lnschema_core.types import FieldAttr

FEATURE_TYPES = {
    "number": "number",
    "int": "int",
    "float": "float",
    "bool": "bool",
    "str": "cat",
    "object": "cat",
}


def convert_numpy_dtype_to_lamin_feature_type(dtype, str_as_cat: bool = True) -> str:
    orig_type = dtype.name
    # strip precision qualifiers
    type = "".join(i for i in orig_type if not i.isdigit())
    if type == "object" or type == "str":
        type = "cat" if str_as_cat else "str"
    return type


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
            if not isinstance(dtype, list) and dtype.__name__ in FEATURE_TYPES:
                dtype_str = FEATURE_TYPES[dtype.__name__]
            else:
                if not isinstance(dtype, list):
                    raise ValueError("dtype has to be a list of Record types")
                registries_str = ""
                for cls in dtype:
                    if not hasattr(cls, "__get_name_with_schema__"):
                        raise ValueError("each element of the list has to be a Record")
                    registries_str += cls.__get_name_with_schema__() + "|"
                dtype_str = f'cat[{registries_str.rstrip("|")}]'
        else:
            dtype_str = dtype
            # add validation that a registry actually exists
            if dtype_str not in FEATURE_TYPES.values() and not dtype_str.startswith(
                "cat"
            ):
                raise ValueError(
                    f"dtype is {dtype_str} but has to be one of 'number', 'int', 'float', 'cat', 'bool', 'cat[...]'!"
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
            categoricals[key] = c
    return categoricals


@classmethod  # type:ignore
@doc_args(Feature.from_df.__doc__)
def from_df(cls, df: pd.DataFrame, field: FieldAttr | None = None) -> RecordsList:
    """{}"""  # noqa: D415
    field = Feature.name if field is None else field
    categoricals = categoricals_from_df(df)

    dtypes = {}
    # categoricals_with_unmapped_categories = {}  # type: ignore
    for name, col in df.items():
        if name in categoricals:
            dtypes[name] = "cat"
        else:
            dtypes[name] = convert_numpy_dtype_to_lamin_feature_type(col.dtype)

    # silence the warning "loaded record with exact same name "
    verbosity = settings.verbosity
    try:
        settings.verbosity = "error"

        registry = field.field.model
        if registry != Feature:
            raise ValueError("field must be a Feature FieldAttr!")
        # create records for all features including non-validated
        features = [Feature(name=name, dtype=dtype) for name, dtype in dtypes.items()]
    finally:
        settings.verbosity = verbosity

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
