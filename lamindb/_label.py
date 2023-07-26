from typing import List, Optional, Union

import pandas as pd
from lamin_utils import logger
from lamindb_setup.dev._docs import doc_args
from lnschema_core import Feature, Label
from lnschema_core.types import ListLike

from lamindb.dev.utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records, index_iterable


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(Label, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 0:
        raise ValueError("Only one non-keyword arg allowed")
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    description: Optional[str] = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    feature: Optional[str] = kwargs.pop("feature") if "feature" in kwargs else None
    feature_id: Optional[str] = (
        kwargs.pop("feature_id") if "feature_id" in kwargs else None
    )
    if len(kwargs) > 0:
        raise ValueError("Only name, description, feature are valid keyword arguments")
    # continue
    if feature is None and feature_id is None:
        logger.warning("Consider passing a corresponding feature for your label!")
    if isinstance(feature, str):
        feature = Feature.select(name=feature).one_or_none()
        if feature is None:
            raise ValueError(
                f"Feature with name {feature} does not exist, please create it:"
                f" ln.Feature(name={feature}, type='float')"
            )
        else:
            feature_id = feature.id
    super(Label, self).__init__(
        name=name, description=description, feature_id=feature_id
    )


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
    "__init__",
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
