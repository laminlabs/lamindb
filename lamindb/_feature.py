from typing import List

from lamindb_setup.dev._docs import doc_args
from lnschema_core import Feature, FeatureValue
from pandas.api.types import is_categorical_dtype, is_string_dtype

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._save import bulk_create


@classmethod  # type:ignore
@doc_args(Feature.from_df.__doc__)
def from_df(cls, df) -> List["Feature"]:
    """{}"""
    records = Feature.from_values(df.columns, field=Feature.name)
    string_or_categorical_columns = [
        col
        for col in df.columns
        if is_string_dtype(df[col]) or is_categorical_dtype(df[col])
    ]
    assert len(records) == len(df.columns)

    for record in records:
        if record.name in string_or_categorical_columns:
            record.type = "str"
            feature = Feature.select(name=record.name).one_or_none()
            values = df[record.name].unique()
            if feature is not None:
                record._values_records = FeatureValue.from_values(
                    values, feature=feature
                )
            else:
                record._values_raw = values
        else:
            record.type = df[record.name].dtype.name
    return records


@doc_args(Feature.save.__doc__)
def save(self, *args, **kwargs) -> None:
    """{}"""
    super(Feature, self).save(*args, **kwargs)
    records = None
    if hasattr(self, "_values_records"):
        records = self._values_records
    if hasattr(self, "_values_raw"):
        records = FeatureValue.from_values(self._values_raw, feature=self)
    if records is not None:
        bulk_create(records)


METHOD_NAMES = [
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
