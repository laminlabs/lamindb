from __future__ import annotations

import lamindb_setup as ln_setup
from lamin_utils import logger

from lamindb.models import ULabel

from ._utils import attach_func_to_class_method


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(ULabel, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 0:
        raise ValueError("Only one non-keyword arg allowed")
    name: str = kwargs.pop("name") if "name" in kwargs else None
    type: str | None = kwargs.pop("type") if "type" in kwargs else None
    is_type: str | None = kwargs.pop("is_type") if "is_type" in kwargs else None
    description: str | None = (
        kwargs.pop("description") if "description" in kwargs else None
    )
    reference: str | None = kwargs.pop("reference") if "reference" in kwargs else None
    reference_type: str | None = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    if len(kwargs) > 0:
        raise ValueError(
            "Only name, description, reference, reference_type are valid keyword arguments"
        )
    if is_type:
        if name.endswith("s"):
            logger.warning(
                "`name` ends with 's', in case you're naming with plural, consider the singular for a type name"
            )
        if name[0].islower():
            logger.warning(
                "`name` starts with lowercase, in case you're naming a type, consider starting with uppercase"
            )
    super(ULabel, self).__init__(
        name=name,
        type=type,
        is_type=is_type,
        description=description,
        reference=reference,
        reference_type=reference_type,
    )


METHOD_NAMES = [
    "__init__",
]

if ln_setup._TESTING:
    from inspect import signature

    SIGS = {
        name: signature(getattr(ULabel, name))
        for name in METHOD_NAMES
        if name != "__init__"
    }

for name in METHOD_NAMES:
    attach_func_to_class_method(name, ULabel, globals())
