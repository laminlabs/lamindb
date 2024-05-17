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
    "int": "int",
    "float": "float",
    "bool": "bool",
    "str": "cat",
    "object": "cat",
}


def convert_numpy_dtype_to_lamin_feature_type(dtype) -> str:
    orig_type = dtype.name
    # strip precision qualifiers
    type = "".join(i for i in orig_type if not i.isdigit())
    if type == "object" or type == "str":
        type = "cat"
    return type


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(Feature, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) != 0:
        raise ValueError("Only non-keyword args allowed")
    type: type | str = kwargs.pop("type") if "type" in kwargs else None
    # cast type
    if type is None:
        raise ValueError("Please pass a type!")
    elif type is not None:
        if not isinstance(type, str):
            if not isinstance(type, list) and type.__name__ in FEATURE_TYPES:
                type_str = FEATURE_TYPES[type.__name__]
            else:
                if not isinstance(type, list):
                    raise ValueError("type has to be a list of Registry types")
                registries_str = ""
                for cls in type:
                    if not hasattr(cls, "__get_name_with_schema__"):
                        raise ValueError(
                            "each element of the list has to be a Registry"
                        )
                    registries_str += cls.__get_name_with_schema__() + "|"
                type_str = f'cat[{registries_str.rstrip("|")}]'
        else:
            type_str = type
            # add validation that a registry actually exists
            if type_str not in FEATURE_TYPES.values() and not type_str.startswith(
                "cat"
            ):
                raise ValueError(
                    "type has to be one of 'number', 'cat', 'bool', 'cat[...]'!"
                )
            if type_str != "cat" and type_str.startswith("cat"):
                registries_str = type_str.replace("cat[", "").rstrip("]")
                if registries_str != "":
                    registry_str_list = registries_str.split("|")
                    for registry_str in registry_str_list:
                        if registry_str not in dict_schema_name_to_model_name(Artifact):
                            raise ValueError(
                                f"'{registry_str}' is an invalid type, pass, e.g. `[ln.ULabel, bt.CellType]` or similar"
                            )
    kwargs["type"] = type_str
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
    """{}."""
    field = Feature.name if field is None else field
    categoricals = categoricals_from_df(df)

    types = {}
    # categoricals_with_unmapped_categories = {}  # type: ignore
    for name, col in df.items():
        if name in categoricals:
            types[name] = "cat"
            # below is a harder feature to write, now, because it requires to
            # query the link tables between the label Registry and file or collection
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

    # silence the warning "loaded record with exact same name "
    verbosity = settings.verbosity
    try:
        settings.verbosity = "error"

        registry = field.field.model
        if registry != Feature:
            raise ValueError("field must be a Feature FieldAttr!")
        # create records for all features including non-validated
        features = [Feature(name=name, type=type) for name, type in types.items()]
    finally:
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
    return RecordsList(features)


# def from_df(
#     self,
#     df: "pd.DataFrame",
#     field: Optional[FieldAttr] = Feature.name,
#     **kwargs,
# ) -> Dict:
#     feature_set = FeatureSet.from_df(df, field=field, **kwargs)
#     if feature_set is not None:
#         feature_sets = {"columns": feature_set}
#     else:
#         feature_sets = {}
#     return feature_sets


@doc_args(Feature.save.__doc__)
def save(self, *args, **kwargs) -> Feature:
    """{}."""
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
