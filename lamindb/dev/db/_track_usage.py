from lndb_setup import settings
from lnschema_core import type

from ._insert import insert


def track_usage(dobject_id, usage_type: type.usage):
    usage_id = insert.usage(  # type: ignore
        type=usage_type,
        user_id=settings.user.id,
        dobject_id=dobject_id,
    )

    return usage_id
