from lnschema_core.ids import base62
from lnschema_core.models import Transform


def __init__(transform, *args, **kwargs):
    if len(args) > 0:  # initialize with all fields from db as args
        super(Transform, transform).__init__(*args, **kwargs)
        return None
    else:  # user-facing calling signature
        # set default ids
        if "id" not in kwargs and "stem_id" not in kwargs:
            kwargs["id"] = base62(14)
            kwargs["stem_id"] = kwargs["id"][:12]
        elif "stem_id" in kwargs:
            assert isinstance(kwargs["stem_id"], str) and len(kwargs["stem_id"]) == 12
            kwargs["id"] = kwargs["stem_id"] + base62(2)
        elif "id" in kwargs:
            assert isinstance(kwargs["id"], str) and len(kwargs["id"]) == 14
            kwargs["stem_id"] = kwargs["id"][:12]
        super(Transform, transform).__init__(**kwargs)


Transform.__init__ = __init__
