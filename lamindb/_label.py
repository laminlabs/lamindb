from typing import List, Optional, Union

import pandas as pd
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Feature, Label
from lnschema_core.types import ListLike

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records, index_iterable


@classmethod  # type:ignore
@doc_args(Label.from_values.__doc__)
def from_values(
    cls, values: ListLike, feature: Optional[Union[Feature, str]] = None, **kwargs
) -> List["Label"]:
    """{}"""
    iterable_idx = index_iterable(values)
    if feature is None and isinstance(values, pd.Series):
        feature = values.name
    if isinstance(feature, str):
        feature = Feature.select(name=feature).one()
    records = get_or_create_records(
        iterable=iterable_idx,
        field=Label.name,
        # here, feature_id is a kwarg, which is an additional condition
        # in queries for potentially existing records
        feature=feature,
    )
    return records


METHOD_NAMES = [
    "from_values",
]

if _TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(Label, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Label, globals())
