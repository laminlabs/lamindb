from __future__ import annotations

from typing import TYPE_CHECKING, Any, get_args

import lamindb_setup as ln_setup
import pandas as pd
from lamin_utils import logger
from lamindb_setup.core._docs import doc_args
from pandas.api.types import CategoricalDtype, is_string_dtype

from lamindb._record import _get_record_kwargs
from lamindb.base.types import FeatureDtype
from lamindb.errors import FieldValidationError, ValidationError
from lamindb.models import Artifact, Feature, Record

from ._query_set import RecordList
from ._utils import attach_func_to_class_method
from .core.relations import dict_module_name_to_model_name

if TYPE_CHECKING:
    from collections.abc import Iterable

    from pandas.core.dtypes.base import ExtensionDtype

    from lamindb.base.types import FieldAttr


FEATURE_DTYPES = set(get_args(FeatureDtype))


def parse_dtype(dtype_str: str) -> list[dict[str, str]]:
    result = []
    # simple dtypes are in FEATURE_DTYPES, composed dtypes are in the form `cat...`
    # if we don't have any of these, throw an error
    if dtype_str not in FEATURE_DTYPES and not dtype_str.startswith("cat"):
        raise ValueError(f"dtype is {dtype_str} but has to be one of {FEATURE_DTYPES}!")
    # now deal with composed categorical dtypes
    if dtype_str != "cat" and dtype_str.startswith("cat"):
        related_registries = dict_module_name_to_model_name(Artifact)
        assert dtype_str.endswith("]")  # noqa: S101
        registries_str = dtype_str.replace("cat[", "")[:-1]  # strip last ]
        if registries_str != "":
            registry_str_list = registries_str.split("|")
            for cat_single_dtype_str in registry_str_list:
                split_result = cat_single_dtype_str.split("[")
                # has sub type
                sub_type_str = ""
                if len(split_result) == 2:
                    registry_str = split_result[0]
                    assert "]" in split_result[1]  # noqa: S101
                    sub_type_field_split = split_result[1].split("].")
                    if len(sub_type_field_split) == 1:
                        sub_type_str = sub_type_field_split[0].strip("]")
                        field_str = ""
                    else:
                        sub_type_str = sub_type_field_split[0]
                        field_str = sub_type_field_split[1]
                elif len(split_result) == 1:
                    registry_field_split = split_result[0].split(".")
                    if (
                        len(registry_field_split) == 2
                        and registry_field_split[1][0].isupper()
                    ) or len(registry_field_split) == 3:
                        # bionty.CellType or bionty.CellType.name
                        registry_str = (
                            f"{registry_field_split[0]}.{registry_field_split[1]}"
                        )
                        field_str = (
                            ""
                            if len(registry_field_split) == 2
                            else registry_field_split[2]
                        )
                    else:
                        # ULabel or ULabel.name
                        registry_str = registry_field_split[0]
                        field_str = (
                            ""
                            if len(registry_field_split) == 1
                            else registry_field_split[1]
                        )
                if registry_str not in related_registries:
                    raise ValidationError(
                        f"'{registry_str}' is an invalid dtype, has to be registry, e.g. ULabel or bionty.CellType"
                    )
                if sub_type_str != "":
                    pass
                    # validate that the subtype is a record in the registry with is_type = True
                registry = related_registries[registry_str]
                if field_str != "":
                    pass
                    # validate that field_str is an actual field of the module
                else:
                    field_str = (
                        registry._name_field
                        if hasattr(registry, "_name_field")
                        else "name"
                    )
                result.append(
                    {
                        "registry": registry,
                        "registry_str": registry_str,
                        "subtype_str": sub_type_str,
                        "field_str": field_str,
                        "field": getattr(registry, field_str),
                    }
                )
    return result


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
            if not hasattr(registry, "__get_name_with_module__"):
                raise ValueError(error_message)
            registries_str += registry.__get_name_with_module__() + "|"
        dtype_str = f"cat[{registries_str.rstrip('|')}]"
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
    name: str = kwargs.pop("name", None)
    dtype: type | str | None = kwargs.pop("dtype") if "dtype" in kwargs else None
    is_type: bool = kwargs.pop("is_type") if "is_type" in kwargs else False
    type_: Feature | str | None = kwargs.pop("type") if "type" in kwargs else None
    if kwargs:
        valid_keywords = ", ".join([val[0] for val in _get_record_kwargs(Feature)])
        raise FieldValidationError(f"Only {valid_keywords} are valid keyword arguments")
    kwargs["name"] = name
    kwargs["type"] = type_
    if is_type:
        if name.endswith("s"):
            logger.warning(
                "`name` ends with 's', in case you're naming with plural, consider the singular for a type name"
            )
        if name[0].islower():
            raise ValidationError(
                "`name` starts with lowercase, name your types with upper-case letters"
            )
        kwargs["is_type"] = is_type
    # cast dtype
    if dtype is None and not is_type:
        raise ValidationError(
            f"Please pass dtype, one of {FEATURE_DTYPES} or a composed categorical dtype"
        )
    dtype_str = None
    if dtype is not None:
        if not isinstance(dtype, str):
            dtype_str = get_dtype_str_from_dtype(dtype)
        else:
            dtype_str = dtype
            parse_dtype(dtype_str)
        kwargs["dtype"] = dtype_str
    super(Feature, self).__init__(*args, **kwargs)
    if not self._state.adding:
        if not (
            self.dtype.startswith("cat") if dtype == "cat" else self.dtype == dtype_str
        ):
            raise ValidationError(
                f"Feature {self.name} already exists with dtype {self.dtype}, you passed {dtype_str}"
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
    registry = field.field.model  # type: ignore
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
        features = [Feature(name=name, dtype=dtype) for name, dtype in dtypes.items()]  # type: ignore
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
