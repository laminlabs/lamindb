from typing import Optional

from lnschema_core.models import Run, Transform


def __init__(run: Run, *args, **kwargs):
    if len(args) == len(run._meta.concrete_fields):
        super(Run, run).__init__(*args, **kwargs)
        return None
    # now we proceed with the user-facing constructor
    if len(args) > 1:
        raise ValueError("Only one non-keyword arg allowed: transform")
    transform: Transform = None
    if "transform" in kwargs or len(args) == 1:
        transform = kwargs.pop("transform") if len(args) == 0 else args[0]
    reference: Optional[str] = (
        kwargs.pop("reference") if "reference" in kwargs else None
    )
    reference_type: Optional[str] = (
        kwargs.pop("reference_type") if "reference_type" in kwargs else None
    )
    if transform is None:
        raise TypeError("Pass transform parameter")
    if transform._state.adding:
        raise ValueError("Please save transform record before creating a run")
    super(Run, run).__init__(
        transform=transform,
        reference=reference,
        reference_type=reference_type,
    )


Run.__init__ = __init__
