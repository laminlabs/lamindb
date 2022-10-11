from lndb_setup import settings
from lnschema_core import type


def track_usage(dobject_id, usage_type: type.usage):
    from ._insert import insert

    usage_id = insert.usage(  # type: ignore
        type=usage_type,
        user_id=settings.user.id,
        dobject_id=dobject_id,
    )

    return usage_id
