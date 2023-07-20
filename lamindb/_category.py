from typing import List

from lamindb_setup.dev._docs import doc_args
from lnschema_core import Category, Feature
from lnschema_core.types import ListLike

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records, index_iterable


@classmethod  # type:ignore
@doc_args(Category.from_values.__doc__)
def from_values(cls, values: ListLike, feature: Feature, **kwargs) -> List["Category"]:
    """{}"""
    iterable_idx = index_iterable(values)
    records = get_or_create_records(
        iterable=iterable_idx,
        field=Category.name,
        # here, feature_id is a kwarg, which is an additional condition
        # in queries for potentially existing records
        feature_id=feature.id,
    )
    return records


METHOD_NAMES = [
    "from_values",
]

if _TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(Category, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, Category, globals())
