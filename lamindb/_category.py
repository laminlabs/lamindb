from typing import List, Union

from lamindb_setup.dev._docs import doc_args
from lnschema_core import Category, Feature
from lnschema_core.types import ListLike

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records, index_iterable


@classmethod  # type:ignore
@doc_args(Category.from_values.__doc__)
def from_values(
    cls, values: ListLike, feature: Union[Feature, str], **kwargs
) -> List["Category"]:
    """{}"""
    iterable_idx = index_iterable(values)
    if isinstance(feature, str):
        feature_id = Feature.select(name=feature).one().id
    else:
        feature_id = feature.id
    records = get_or_create_records(
        iterable=iterable_idx,
        field=Category.name,
        # here, feature_id is a kwarg, which is an additional condition
        # in queries for potentially existing records
        feature_id=feature_id,
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
