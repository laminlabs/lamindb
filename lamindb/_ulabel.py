from __future__ import annotations

import lamindb_setup as ln_setup

from lamindb.models import ULabel

from ._utils import attach_func_to_class_method


def __init__(self, *args, **kwargs):
    if len(args) == len(self._meta.concrete_fields):
        super(ULabel, self).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 0:
        raise ValueError("Only one non-keyword arg allowed")
    name: str | None = kwargs.pop("name") if "name" in kwargs else None
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
    super(ULabel, self).__init__(
        name=name,
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
