import hashlib

from lnschema_core.ids import base62
from lnschema_core.models import Transform

from .dev.hashing import to_b64_str


def __init__(transform, *args, **kwargs):
    if len(args) > 0:  # initialize with all fields from db as args
        super(Transform, transform).__init__(*args, **kwargs)
        return None
    else:  # user-facing calling signature
        if "version" not in kwargs:
            kwargs["version"] = "0"
        elif not isinstance(kwargs["version"], str):
            raise ValueError("version must be str, e.g., '0', '1', etc.")
        id_ext = to_b64_str(hashlib.md5(kwargs["version"].encode()).digest())[:2]
        # set default ids
        if "id" not in kwargs and "stem_id" not in kwargs:
            kwargs["stem_id"] = base62(12)
            kwargs["id"] = kwargs["stem_id"] + id_ext
        elif "stem_id" in kwargs:
            assert isinstance(kwargs["stem_id"], str) and len(kwargs["stem_id"]) == 12
            kwargs["id"] = kwargs["stem_id"] + id_ext
        elif "id" in kwargs:
            assert isinstance(kwargs["id"], str) and len(kwargs["id"]) == 14
            kwargs["stem_id"] = kwargs["id"][:12]
        super(Transform, transform).__init__(**kwargs)


Transform.__init__ = __init__
