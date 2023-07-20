from typing import List

import pandas as pd
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Category, Feature
from pandas.api.types import is_categorical_dtype, is_string_dtype

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._save import bulk_create


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
    records = Feature.from_values(df.columns, field=Feature.name)
    assert len(records) == len(df.columns)

    string_cols = [col for col in df.columns if is_string_dtype(df[col])]
    categoricals = {col: df[col] for col in df.columns if is_categorical_dtype(df[col])}
    for key in string_cols:
        c = pd.Categorical(df[key])
        # TODO: We should only check if non-null values are unique, but
        # this would break cases where string columns with nulls could
        # be written as categorical, but not as string.
        # Possible solution: https://github.com/scverse/anndata/issues/504
        if len(c.categories) < len(c):
            categoricals[key] = c

    for record in records:
        if record.name in categoricals:
            record.type = "Category"
            feature = Feature.select(name=record.name).one_or_none()
            categories = categoricals[record.name].categories
            if feature is not None:
                record._categories_records = Category.from_values(
                    categories, feature=feature
                )
            else:
                record._categories_raw = categories
        else:
            orig_type = df[record.name].dtype.name
            # strip precision qualifiers
            record.type = "".join(i for i in orig_type if not i.isdigit())
    return records


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
