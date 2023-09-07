from typing import Dict, List, Optional, Union

import pandas as pd
from lamindb_setup.dev._docs import doc_args
from lnschema_core.models import Feature, Registry
from pandas.api.types import CategoricalDtype, is_string_dtype

from lamindb._utils import attach_func_to_class_method
from lamindb.dev._settings import settings

from . import _TESTING
from .dev._priors import NUMBER_TYPE


def convert_numpy_dtype_to_lamin_feature_type(dtype) -> str:
    orig_type = dtype.name
    # strip precision qualifiers
    type = "".join(i for i in orig_type if not i.isdigit())
    if type == "int" or type == "float":
        type = NUMBER_TYPE
    return type


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(Feature, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) != 0:
        raise ValueError("Only non-keyword args allowed")
    type: Optional[Union[type, str]] = kwargs.pop("type") if "type" in kwargs else None
    registries: Optional[List[Registry]] = (
        kwargs.pop("registries") if "registries" in kwargs else None
    )
    # cast type
    type_str = None
    if type is not None:
        type_str = type.__name__ if not isinstance(type, str) else type
    kwargs["type"] = type_str
    # cast registries
    registries_str: Optional[str] = None
    if registries is not None:
        if isinstance(registries, str):
            # TODO: add more validation
            registries_str = registries
        else:
            if not isinstance(registries, List):
                raise ValueError("registries has to be a list of Registry types")
            registries_str = ""
            for cls in registries:
                if not hasattr(cls, "__get_name_with_schema__"):
                    raise ValueError("each element of the list has to be a Registry")
                registries_str += cls.__get_name_with_schema__() + "|"
            registries_str = registries_str.rstrip("|")
    kwargs["registries"] = registries_str
    super(Feature, self).__init__(*args, **kwargs)


def categoricals_from_df(df: "pd.DataFrame") -> Dict:
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
def from_df(cls, df: "pd.DataFrame") -> List["Feature"]:
    """{}"""
    categoricals = categoricals_from_df(df)

    types = {}
    # categoricals_with_unmapped_categories = {}  # type: ignore
    for name, col in df.items():
        if name in categoricals:
            types[name] = "category"
            # below is a harder feature to write, now, because it requires to
            # query the link tables between the label Registry and file or dataset
            # the original implementation fell short
            # categorical = categoricals[name]
            # if hasattr(
            #     categorical, "cat"
            # ):  # because .categories > pd2.0, .cat.categories < pd2.0
            #     categorical = categorical.cat
            # categories = categorical.categories
            # categoricals_with_unmapped_categories[name] = ULabel.filter(
            #     feature=name
            # ).inspect(categories, "name", logging=False)["not_mapped"]
        else:
            types[name] = convert_numpy_dtype_to_lamin_feature_type(col.dtype)

    # silence the info "loaded record with exact same name "
    verbosity = settings.verbosity
    settings.verbosity = 1
    # create records for all features including non-validated
    features = [Feature(name=name, type=type) for name, type in types.items()]
    settings.verbosity = verbosity

    assert len(features) == len(df.columns)

    # if len(categoricals_with_unmapped_categories) > 0:
    #     n_max = 20
    #     categoricals_with_unmapped_categories_formatted = "\n      ".join(
    #         [
    #             (
    #                 f"{key} ({len(value)}): {', '.join(value)}"
    #                 if len(value) <= 5
    #                 else f"{key} ({len(value)}): {', '.join(value[:5])} ..."
    #             )
    #             for key, value in take(
    #                 n_max, categoricals_with_unmapped_categories.items()
    #             )
    #         ]
    #     )
    #     if len(categoricals_with_unmapped_categories) > n_max:
    #         categoricals_with_unmapped_categories_formatted += "\n      ..."
    #     categoricals_with_unmapped_categories_formatted
    #     logger.info(
    #         f"{len(categoricals_with_unmapped_categories)} features have"
    #         f" {colors.yellow('unmapped categories')}:\n     "
    #         f" {categoricals_with_unmapped_categories_formatted}"
    #     )
    return features


@doc_args(Feature.save.__doc__)
def save(self, *args, **kwargs) -> None:
    """{}"""
    super(Feature, self).save(*args, **kwargs)


METHOD_NAMES = [
    "__init__",
    "from_df",
    "save",
]

if _TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(Feature, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Feature, globals())
