from typing import List, Optional

from lamindb_setup.dev._docs import doc_args
from lnschema_core import ULabel
from lnschema_core.types import ListLike

from lamindb._utils import attach_func_to_class_method

from . import _TESTING
from ._from_values import get_or_create_records


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(ULabel, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 0:
        raise ValueError("Only one non-keyword arg allowed")
    name: Optional[str] = kwargs.pop("name") if "name" in kwargs else None
    description: Optional[str] = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    reference: Optional[str] = (
        kwargs.pop("reference") if "reference" in kwargs else None
    )
    reference_type: Optional[str] = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    if len(kwargs) > 0:
        raise ValueError(
            "Only name, description, reference, reference_type are valid keyword arguments"  # noqa
        )
    super(ULabel, self).__init__(
        name=name,
        description=description,
        reference=reference,
        reference_type=reference_type,
    )


@classmethod  # type:ignore
@doc_args(ULabel.from_values.__doc__)
def from_values(cls, values: ListLike, **kwargs) -> List["ULabel"]:
    """{}"""
    records = get_or_create_records(
        iterable=values,
        field=ULabel.name,
    )
    return records


METHOD_NAMES = [
    "__init__",
    "from_values",
]

if _TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(ULabel, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, ULabel, globals())
