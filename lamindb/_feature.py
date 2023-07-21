from itertools import islice
from typing import List

import pandas as pd
from lamin_utils import logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Category, Feature
from pandas.api.types import is_categorical_dtype, is_string_dtype

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._save import bulk_create


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
def from_df(cls, df) -> List["Feature"]:
    """{}"""
    features = Feature.from_values(df.columns, field=Feature.name)
    assert len(features) == len(df.columns)

    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {col: df[col] for col in df.columns if is_categorical_dtype(df[col])}
    for key in string_cols:
        c = pd.Categorical(df[key])
        if len(c.categories) < len(c):
            categoricals[key] = c

    categoricals_with_unmapped_categories = {}
    for feature in features:
        if feature.name in categoricals:
            feature.type = "Category"
            categorical = categoricals[feature.name]
            if hasattr(
                categorical, "cat"
            ):  # because .categories > pd2.0, .cat.categories < pd2.0
                categorical = categorical.cat
            categories = categorical.categories
            categoricals_with_unmapped_categories[feature.name] = Category.select(
                feature=feature
            ).inspect(categories, "name", logging=False)["not_mapped"]
        else:
            orig_type = df[feature.name].dtype.name
            # strip precision qualifiers
            feature.type = "".join(i for i in orig_type if not i.isdigit())
    if len(categoricals) > 0:
        categoricals_with_unmapped_categories_formatted = "\n      ".join(
            [
                f"{key}: {', '.join(value)}"
                for key, value in take(7, categoricals_with_unmapped_categories.items())
            ]
        )
        if len(categoricals_with_unmapped_categories) > 7:
            categoricals_with_unmapped_categories_formatted += "\n      ..."
        categoricals_with_unmapped_categories_formatted
        logger.info(
            "There are unmapped categories:\n     "
            f" {categoricals_with_unmapped_categories_formatted}"
        )
        hint_formatted = "\n      ".join(
            [
                f"ln.Category.from_values(df['{key}'])"
                for key in take(7, categoricals_with_unmapped_categories)
            ]
        )
        if len(categoricals_with_unmapped_categories) > 7:
            hint_formatted += "\n      ..."
        logger.hint(f"Consider adding them via:\n      {hint_formatted}")
    return features


@doc_args(Feature.save.__doc__)
def save(self, *args, **kwargs) -> None:
    """{}"""
    super(Feature, self).save(*args, **kwargs)
    records = None
    if hasattr(self, "_categories_records"):
        records = self._categories_records
    if hasattr(self, "_categories_raw"):
        records = Category.from_values(self._categories_raw, feature=self)
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
