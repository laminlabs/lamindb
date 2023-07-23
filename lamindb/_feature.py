from itertools import islice
from typing import List

import pandas as pd
from lamin_utils import logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Feature, Label
from pandas.api.types import is_categorical_dtype, is_string_dtype

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._save import bulk_create


def convert_numpy_dtype_to_lamin_feature_type(dtype) -> str:
    orig_type = dtype.name
    # strip precision qualifiers
    type = "".join(i for i in orig_type if not i.isdigit())
    return type


def take(n, iterable):
    """Return the first n items of the iterable as a list."""
    return list(islice(iterable, n))


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(Feature, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) != 0:
        raise ValueError("Only non-keyword args allowed")
    super(Feature, self).__init__(*args, **kwargs)


@classmethod  # type:ignore
@doc_args(Feature.from_df.__doc__)
def from_df(cls, df: "pd.DataFrame") -> List["Feature"]:
    """{}"""
    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {col: df[col] for col in df.columns if is_categorical_dtype(df[col])}
    for key in string_cols:
        c = pd.Categorical(df[key])
        if len(c.categories) < len(c):
            categoricals[key] = c

    types = {}
    categoricals_with_unmapped_categories = {}
    for name, col in df.items():
        if name in categoricals:
            types[name] = "categorical"
            categorical = categoricals[name]
            if hasattr(
                categorical, "cat"
            ):  # because .categories > pd2.0, .cat.categories < pd2.0
                categorical = categorical.cat
            categories = categorical.categories
            categoricals_with_unmapped_categories[name] = Label.select(
                feature=name
            ).inspect(categories, "name", logging=False)["not_mapped"]
        else:
            types[name] = convert_numpy_dtype_to_lamin_feature_type(col.dtype)

    features = Feature.from_values(df.columns, field=Feature.name, types=types)
    assert len(features) == len(df.columns)

    if len(categoricals) > 0:
        n_max = 20
        categoricals_with_unmapped_categories_formatted = "\n      ".join(
            [
                f"{key}: {', '.join(value)}"
                for key, value in take(
                    n_max, categoricals_with_unmapped_categories.items()
                )
            ]
        )
        if len(categoricals_with_unmapped_categories) > n_max:
            categoricals_with_unmapped_categories_formatted += "\n      ..."
        categoricals_with_unmapped_categories_formatted
        logger.info(
            "There are unmapped categories:\n     "
            f" {categoricals_with_unmapped_categories_formatted}"
        )
    return features


@doc_args(Feature.save.__doc__)
def save(self, *args, **kwargs) -> None:
    """{}"""
    super(Feature, self).save(*args, **kwargs)
    records = None
    if hasattr(self, "_categories_records"):
        records = self._categories_records
    if hasattr(self, "_categories_raw"):
        records = Label.from_values(self._categories_raw, feature=self)
    if records is not None:
        bulk_create(records)


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
